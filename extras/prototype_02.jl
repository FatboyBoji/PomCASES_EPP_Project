using JuMP
using HiGHS

"""
### Hydrogen Management Optimization Model
"""

function SolveEPP(time_limit::Int64)
    #######################
    ### Model Definition ###
    #######################

    EPP = Model(HiGHS.Optimizer)
    set_optimizer_attribute(EPP, "time_limit", time_limit)

    #######################
    ### Data Parameters ###
    #######################

    # Existing data
    J = 15; #Number of jobs
    M = 2; #Number of machines
    P = 24; #Number of periods
    S = 5; #Number of states; # 1:off 2:turnon 3:idle 4:process 5:turnoff

    # Processing time of job j on machine m
    pt = [1 1 2 3 3 1 1 1 2 3 2 3 3 1 1
          2 2 4 6 6 2 2 2 4 6 2 4 6 6 2]; 
          
    # Ramp up time of machine m
    ps = [2 2]; 

    # Maximal setup time per machine
    BigK = ps; 
      
    # Index sets of possible state transitions
    T = [[1,2],     # successors of state 1 (off)
        [2,3,4],    # successors of state 2 (turn on)
        [3,4,5],    # successors of state 3 (idle)
        [3,4,5],    # successors of state 4 (processing)
        [5,1]       # successors of state 5 (turn off)
        ]; 

    # Energy demand per machine and state
    e = [0 190 10 100 20  # machine 1
         0  90 10  50 10   ];  # machine 2

    # Peak load limit 
    c_PL = 200;

    # Energy prices per period
    c_buy = [1, 1, 1, 1, 2, 2, 3, 4, 5, 5, 1, 4, 3, 3, 1, 1, 1, 1, 4, 3, 3, 1, 1, 1];
    c_sell = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];

    # Renewable generation per period
    re = [0, 0, 0, 20, 30, 40, 40, 50, 70, 80, 90, 80, 100, 70, 70, 60, 50, 30, 30, 10, 10, 0, 0, 0];

    # Battery-related parameters (inactive by default)
    Max_charge_rate = 0;     
    Max_discharge_rate = 0;  
    Bat_cap = 0; 
    Bat_SoC_initial = 0; 
    Bat_SoC_end = 0; 

    # Add cost parameters
    c_charge = 1.0;    # €/kWh cost of charging
    c_discharge = 0.5; # €/kWh cost of discharging

    # Battery efficiency parameters
    eta_charge = 0;    # Charging efficiency
    eta_discharge = 0; # Discharging efficiency
    Bat_autodischarge_rate = 0;  # Self-discharge rate per period

    ###############  Decision variables   ###############
    #Scheduling
    @variable(EPP, X[1:M,1:J,1:P], Bin); #machine m processes job j in period p
    @variable(EPP, Y[1:M,1:J,1:P], Bin); #machine m begins job j in period p
    @variable(EPP, alpha[1:M,1:P,1:S], Bin); #machine m is in period p in state s

    #Energy
    @variable(EPP, E_sell[1:P] >= 0); #Energy sold in kWh in period p
    @variable(EPP, E_buy[1:P] >= 0); #Energy purchased in kWh in period p
    
    @variable(EPP, E_charge[1:P] >= 0); #Energy send to battery in kWh in period p 
    @variable(EPP, E_discharge[1:P] >= 0); #Energy taken from battery in kWh in period p 
    @variable(EPP, Bat_SoC[1:P] >= 0); #Battery state of charge in period p 

    #######################
    ### old constraints ...
    #######################
    # (2) machine must be in operating mode (s=4) when a job is produced 
    @constraint(EPP, [m=1:M, p=1:P], 
    	sum(X[m,j,p] for j=1:J) == alpha[m,p,4]); 

    # (3) machine must be in excatly one state at a time
    @constraint(EPP, [m=1:M,p=1:P], 
        sum(alpha[m,p,s] for s=1:S) == 1); 

    # (4) only defined state transistions allowed
    @constraint(EPP, [m=1:M,p=1:P-1,s=1:S], 
    alpha[m,p,s] <= sum(alpha[m,p+1,s_prime] for s_prime=T[s])); 

    # (5) no preemption: processing of jobs cannot be interrupted
    @constraint(EPP, [m=1:M,j=1:J, p=1:(P-pt[m,j]+1)], 
        sum(X[m,j,p_prime] for p_prime=p:(p+pt[m,j]-1)) >= pt[m,j]*Y[m,j,p]); 

    # (6) Setup constraint: production only after initial setup of the machine 
    @constraint(EPP, [m=1:M,p=ps[m]+1:P], 
    (alpha[m,p-1,3]+alpha[m,p-1,4])*BigK[m] + sum(alpha[m,p_prime,2] for p_prime=p-ps[m]:p-1) >= ps[m]*sum(Y[m,j,p] for j=1:J)); 

    # (7) Uninterrupted setup
    @constraint(EPP, [m=1:M, p=ps[m]:P-1],
        sum(alpha[m,p_prime,2] for p_prime=(p-ps[m]+1):p) >= BigK[m]*(alpha[m,p,2]-alpha[m,p+1,2]));

    # (8) every job starts once on one machine
    @constraint(EPP, [j=1:J], 
        sum(Y[m,j,p] for m=1:M, p=(ps[m]+1):(P-pt[m,j]+1)) == 1); 

    # (9) each machine can process only one job at a time
    @constraint(EPP, [m=1:M,p=ps[m]:P], 
        sum(X[m,j,p] for j=1:J) <= 1); 

    # (10) machine is turned off or in setup mode for the first periods
    @constraint(EPP, [m=1:M, p=1:ps[m]], 
        alpha[m,p,1]+alpha[m,p,2] == 1); 

    # (11) machine is turned off in last period
    @constraint(EPP, [m=1:M, p=P:P], 
        alpha[m,p,1] == 1); 

    # (12) energy balance
    @constraint(EPP, [p=1:P], 
        sum(alpha[m,p,s]*e[m,s] for m=1:M, s=1:S) + E_sell[p] + E_charge[p] == E_buy[p] + re[p] + E_discharge[p]);  

    # (13) Peak load
    @constraint(EPP, [p=1:P], 
        c_PL>=E_sell[p]+E_buy[p]);

    # (14), (15) charge and discharge constraints 
    @constraint(EPP, [p=1:P], 
        E_charge[p]<=Max_charge_rate);

    @constraint(EPP, [p=1:P], 
        E_discharge[p]<=Max_discharge_rate);

    # (16) Battery state of charge 
    @constraint(EPP, [p=1:1], # first period
        Bat_SoC[p]==Bat_SoC_initial+E_charge[p]-E_discharge[p]);

    @constraint(EPP, [p=2:P], # all other periods
        Bat_SoC[p]==Bat_SoC[p-1]+E_charge[p]-E_discharge[p]);

    @constraint(EPP, [p=P:P], # last period
        Bat_SoC_end==Bat_SoC[p-1]+E_charge[p]-E_discharge[p]);

    # (17) Battery capacity 
    @constraint(EPP, [p=1:P], 
            Bat_SoC[p]<=Bat_cap);


    ############################################
    ### Hydrogen fuel cell 
    ############################################
    H_max = 0; 
    H_storage_initial = 0; 
    H_storage_end = 0;      
    H_Max_charge_rate = 0;     
    H_Max_discharge_rate = 0;

    #euro/kg wasserstoffpreise
    c_hBuy = [1, 1, 1, 1, 2, 2, 3, 4, 5, 5, 1, 4, 3, 3, 1, 1, 1, 1, 4, 3, 3, 1, 1, 1]
    c_hSell = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    #lagerkosten 20-750 $
    c_H_storage = 400

    H_heizwert = 1/3; # kg/kwh
    
    a_toH2 = 0.63
    a_fromH2 = 0.48
    
    leistung_punkte              = [0,   10,   20,   40,   60,   80,  100]                 # kWh (Leistung)
    wirkungsgrad_punkte_toH2     = [0, 0.33, 0.45, 0.55, 0.58, 0.52, 0.45]    # toH2
    wirkungsgrad_punkte_fromH2   = [0, 0.33, 0.45, 0.55, 0.58, 0.52, 0.45]    # fromH2


    #######################
    ### Decision Variables
    #######################          
    E_toH2 >= 0
    E_fromH2 >= 0
    H_storage >= 0
    H_buy[p] >= 0   #kg
    H_cell[p] >= 0  #kg


    # y_elec[p] ∈ {0,1}    # Elektrolyse-Modus aktiv
    # y_fc[p] ∈ {0,1}      # Brennstoffzellen-Modus aktiv
    # y_idle[p] ∈ {0,1}    # Leerlauf-Modus aktiv

    #sos2 wirkungsgrade
    # SOS2-basierte Wirkungsgrade:
    # 1. Für Elektrolyse (toH2)
    E_toH2[p] = Σ(λ_toH2[p,i] * leistung_punkte[i])     # Leistungsinterpolation
    η_toH2[p] = Σ(λ_toH2[p,i] * wirkungsgrad_punkte_toH2[i])   # Wirkungsgradinterpolation
    Σ(λ_toH2[p,i]) = 1
    λ_toH2[p,:] ∈ SOS2

    # 2. Für Brennstoffzelle (fromH2)
    E_fromH2[p] = Σ(λ_fromH2[p,i] * leistung_punkte[i])   # Leistungsinterpolation
    η_fromH2[p] = Σ(λ_fromH2[p,i] * wirkungsgrad_punkte_fromH2[i])   # Wirkungsgradinterpolation
    Σ(λ_fromH2[p,i]) = 1
    λ_fromH2[p,:] ∈ SOS2

    EC_menge[p] = E_toH2[p] * η_toH2[p] * H_heizwert
    FC_menge[p] = E_fromH2[p] * η_fromH2[p] * H_heizwert

    #######################
    ### Constraints + sos2 ,   EC - Elektrolyse,   FC - Fuelcell
    #######################
    # max constraint:
    H_storage[p] <= H_max
    E_toH2[p] ≤ H_Max_charge_rate
    E_fromH2[p] ≤ H_Max_discharge_rate
    E_toH2[p] * E_fromH2[p] == 0  # Entweder Elektrolyse oder Brennstoffzelle

    # speicher 
    H_storage[1] = H_storage_initial + EC[p] - FC[p] + H_buy[p] - H_sell[p]

    H_storage[p] = H_storage[p-1] + EC_menge[p] - FC_menge[p] + H_buy[p] - H_sell[p]

    H_storage_end[p] = H_storage[p-1] + EC_menge[p] - FC_menge[p] + H_buy[p] - H_sell[p]

    ## 
    EC_menge[p] = E_toH2[p]   * a_toH2   * H_heiwert      # kwh * % * kg/kwh = kg
    FC_menge[p] = E_fromH2[p] * a_fromH2 * H_heizwert


    # !!! new energy balance
    @constraint(EPP, [p=1:P], 
        sum(alpha[m,p,s]*e[m,s] for m=1:M, s=1:S) + E_sell[p] + E_charge[p] + E_toH2 == E_buy[p] + re[p] + E_discharge[p] + E_fromH2);  

    """
    1kg H2 -> 100% -> xkwh
    1kg H2 -> 50% -> ykwh

    40% -> 
    50% -> 50kwH 2 kg
    60% -> 50kwH 2.5kg H2;
    """

    #######################
    ### old without the fuel cell Objective Function 
    #######################
    @objective(EPP, Min, 
        sum(c_buy[p]*E_buy[p] - c_sell[p]*E_sell[p] for p=1:P) + # original energy costs
        sum(c_charge*E_charge[p] + c_discharge*E_discharge[p] for p=1:P), # battery operation costs
        sum(c_hBuy[p] * H_buy[p] + c_hSell[p] + H_sell[p] for p=1:P),
        sum(c_H_storage * H_storage[p] for p=1:P)
    );

    #######################
    ### Solving the Model ###
    #######################

    optimize!(EPP)

    #######################
    ### Return Results ###
    #######################

    return EPP
end

####################
# Main Code Section
####################

KPI = SolveEPP(360);
