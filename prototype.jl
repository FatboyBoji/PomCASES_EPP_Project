using JuMP, Gurobi, Statistics, StatsPlots, DataFrames, CSV, Dates, Printf


"""
### Hydrogen Management Optimization Model
"""

function SolveEPP(time_limit::Int64)
    #######################
    ### Model Definition ###
    #######################

    EPP = Model(Gurobi.Optimizer)
    set_optimizer_attribute(EPP, "TimeLimit", time_limit)
    set_optimizer_attribute(EPP, "OutputFlag", 0)

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
    c_sell = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];

    # Renewable generation per period
    re = [0, 0, 0, 20, 30, 40, 40, 50, 70, 80, 90, 80, 100, 70, 70, 60, 50, 30, 30, 10, 10, 0, 0, 0];

    # Battery-related parameters (inactive by default)
    Max_charge_rate = 0 #40
    Max_discharge_rate = 0 #40
    Bat_cap = 0 #300
    Bat_SoC_initial = 0 #150
    Bat_SoC_end = 0 #150

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
    # @constraint(EPP, [p=1:P], 
    #     sum(alpha[m,p,s]*e[m,s] for m=1:M, s=1:S) + E_sell[p] + E_charge[p] == E_buy[p] + re[p] + E_discharge[p]);  

    # (13) Peak load
    @constraint(EPP, [p=1:P], c_PL>=E_sell[p]+E_buy[p]);

    # (14), (15) charge and discharge constraints 
    @constraint(EPP, [p=1:P], E_charge[p]<=Max_charge_rate);

    @constraint(EPP, [p=1:P], E_discharge[p]<=Max_discharge_rate);

    # (16) Battery state of charge 
    @constraint(EPP, [p=1:1], # first period
        Bat_SoC[p]==Bat_SoC_initial+E_charge[p]-E_discharge[p]);

    @constraint(EPP, [p=2:P], # all other periods
        Bat_SoC[p]==Bat_SoC[p-1]+E_charge[p]-E_discharge[p]);

    @constraint(EPP, [p=P:P], # last period
        Bat_SoC_end==Bat_SoC[p-1]+E_charge[p]-E_discharge[p]);

    # (17) Battery capacity 
    @constraint(EPP, [p=1:P], Bat_SoC[p]<=Bat_cap);


    ############################################
    ### Hydrogen fuel cell 
    ############################################
    H_max = 1000;                # kg H2 Speicherkapazität
    H_storage_initial = 0;       # kg H2 Anfangsbestand
    H_storage_end = 0;           # kg H2 Endbestand     
    H_Max_charge_rate = 100;     # kW max Elektrolyse-Leistung    
    H_Max_discharge_rate = 100;  # kW max Brennstoffzellen-Leistung


    #euro/kg wasserstoffpreise
    c_hBuy = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
    c_hSell = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]

    #lagerkosten 20-750 $ per kg H2
    c_H_storage = 5

    #kosten für die umwandlung
    c_h2charge  = 0.15      # €/kWh cost of charging
    c_h2discharge = 0.15    # €/kWh cost of charging

    H_heizwert = 1/3; # kg/kwh //  33.33 kwh/kg
    
    # a_toH2 = 0.63
    # a_fromH2 = 0.48
    
    leistung_punkte              = [0,   10,    20,    30,    40,    50,   60,    70,   80,  90,  100]   # kW (Leistung)
    wirkungsgrad_punkte_toH2     = [0, 0.75,  0.78,  0.72,  0.65,  0.65,  0.6,   0.6, 0.55, 0.5,  0.5]   # toH2
    wirkungsgrad_punkte_fromH2   = [0,  0.6,  0.58,  0.55,  0.52,  0.50, 0.48,  0.45, 0.43, 0.4,  0.3]   # fromH2

    #######################
    ### Decision Variables ###
    #######################     
    # variable für "wenn storage[p] < max dann EC möglich ,  storage > 0 dann FC möglich"     
    @variable(EPP, E_toH2[1:P] >= 0)
    @variable(EPP, E_fromH2[1:P] >= 0)
    @variable(EPP, H_storage[1:P] >= 0)
    @variable(EPP, H_buy[1:P] >= 0)   #kg
    @variable(EPP, H_sell[1:P] >= 0)  #kg
    
    # SOS2 Variablen
    @variable(EPP, λ_toH2[1:P, 1:length(leistung_punkte)] >= 0)
    @variable(EPP, λ_fromH2[1:P, 1:length(leistung_punkte)] >= 0)
    @variable(EPP, η_toH2[1:P] >= 0)
    @variable(EPP, η_fromH2[1:P] >= 0)
    @variable(EPP, EC_menge[1:P] >= 0)
    @variable(EPP, FC_menge[1:P] >= 0)

    @variable(EPP, active_1[1:P], Bin)
    @variable(EPP, active_2[1:P], Bin)
    @variable(EPP, inv_η_toH2[1:P] >= 0)
    @variable(EPP, inv_η_fromH2[1:P] >= 0)

    #######################
    ### Constraints + sos2 ,   EC - Elektrolyse,   FC - Fuelcell
    #######################
    # max constraint:
    @constraint(EPP, [p=1:P], H_storage[p] <= H_max)
    @constraint(EPP, [p=1:P], E_toH2[p] <= H_Max_charge_rate)
    @constraint(EPP, [p=1:P], E_fromH2[p] <= H_Max_discharge_rate)

    # @constraint(EPP, [p=1:P], inv_η_toH2[p] * η_toH2[p] == 1)
    # @constraint(EPP, [p=1:P], inv_η_fromH2[p] * η_fromH2[p] == 1)

    # SOS2 Constraints für Elektrolyse
    @constraint(EPP, [p=1:P], 
        E_toH2[p] == sum(λ_toH2[p,i] * leistung_punkte[i] for i in eachindex(leistung_punkte)) * active_1[p])
    @constraint(EPP, [p=1:P], active_1[p] <= 1.999 - (H_storage[p] / H_max) )  # Aktivieren nur wenn H_storage < H_max
    # @constraint(EPP, [p=1:P], 
    #     E_toH2[p] == sum(λ_toH2[p,i] * leistung_punkte[i] for i in eachindex(leistung_punkte)))
    @constraint(EPP, [p=1:P], 
        η_toH2[p] == sum(λ_toH2[p,i] * wirkungsgrad_punkte_toH2[i] for i in eachindex(leistung_punkte)))
    @constraint(EPP, [p=1:P], 
        sum(λ_toH2[p,i] for i in eachindex(leistung_punkte)) == 1)
    @constraint(EPP, [p=1:P], λ_toH2[p,:] in SOS2())

    # SOS2 Constraints für Brennstoffzelle
    @constraint(EPP, [p=1:P], E_fromH2[p] == sum(λ_fromH2[p,i] * leistung_punkte[i] for i in eachindex(leistung_punkte)) * active_2[p])
    @constraint(EPP, [p=1:P], active_2[p] <= H_storage[p])  # Aktivieren nur wenn H_storage > 0
    # @constraint(EPP, [p=1:P], 
    #     E_fromH2[p] == sum(λ_fromH2[p,i] * leistung_punkte[i] for i in eachindex(leistung_punkte)))
    @constraint(EPP, [p=1:P], 
        η_fromH2[p] == sum(λ_fromH2[p,i] * wirkungsgrad_punkte_fromH2[i] for i in eachindex(leistung_punkte)))
    @constraint(EPP, [p=1:P], sum(λ_fromH2[p,i] for i in eachindex(leistung_punkte)) == 1)
    @constraint(EPP, [p=1:P], λ_fromH2[p,:] in SOS2())

    # Umrechnung Energie zu H2-Menge
    @constraint(EPP, [p=1:P], 
        EC_menge[p] == E_toH2[p] * η_toH2[p] * H_heizwert)  # kwh * % * kg/kwh = kg
    @constraint(EPP, [p=1:P], 
        FC_menge[p] == E_fromH2[p] * η_fromH2[p] * H_heizwert) 

    ## wenn storage[p] < max dann EC möglich ,  storage > 0 dann FC möglich

    # Speicherstandsbilanz
    # first period
    @constraint(EPP, [p=1:1], 
        H_storage[p] == H_storage_initial + EC_menge[p] - FC_menge[p] + H_buy[p] - H_sell[p])
    
    # other periods
    @constraint(EPP, [p=2:P],
        H_storage[p] == H_storage[p-1] + EC_menge[p] - FC_menge[p] + H_buy[p] - H_sell[p])
    
    # last period
    @constraint(EPP, [p=P:P], 
        H_storage_end == H_storage[p-1] + EC_menge[p] - FC_menge[p] + H_buy[p] - H_sell[p])

    # eine bedingung damit EC und FC nicht gleichzeitig laufen
    @constraint(EPP, [p=1:P], E_toH2[p] * E_fromH2[p] == 0)


    ### kosten timer preis beim umschalten von EC zu FC

    # Energiebilanz
    @constraint(EPP, [p=1:P], 
        sum(alpha[m,p,s]*e[m,s] for m=1:M, s=1:S) + E_sell[p] + E_charge[p] + E_toH2[p] == 
        E_buy[p] + re[p] + E_discharge[p] + E_fromH2[p])

    #######################
    ### Objective Function ###
    #######################
    function safe_inverse(x)
        return x == 0 ? 1000.0 : 1/x  
    end
    
    wirkungsgrad_punkte_toH2_inv = [safe_inverse(x) for x in wirkungsgrad_punkte_toH2]
    wirkungsgrad_punkte_fromH2_inv = [safe_inverse(x) for x in wirkungsgrad_punkte_fromH2]


    @constraint(EPP, [p=1:P], inv_η_toH2[p] == sum(λ_toH2[p,i] * wirkungsgrad_punkte_toH2_inv[i] for i in eachindex(leistung_punkte)))
    @constraint(EPP, [p=1:P], inv_η_fromH2[p] == sum(λ_fromH2[p,i] * wirkungsgrad_punkte_fromH2_inv[i] for i in eachindex(leistung_punkte)))

    @objective(EPP, Min, 
        sum(c_buy[p]*E_buy[p] - c_sell[p]*E_sell[p] for p=1:P) +  
        sum(c_charge*E_charge[p] + c_discharge*E_discharge[p] for p=1:P) +  
        sum(c_h2charge * inv_η_toH2[p] * E_toH2[p] + c_h2discharge * inv_η_fromH2[p] * E_fromH2[p] for p=1:P) +
        sum(c_hBuy[p] * H_buy[p] - c_hSell[p] * H_sell[p] for p=1:P) +  
        sum(c_H_storage * H_storage[p] for p=1:P)  
    )

    #######################
    ### Solving the Model ###
    #######################

    optimize!(EPP)

    #######################
    ### Return Results ###
    #######################

    
    ED = [sum(value.(alpha[m,p,s])*e[m,s] for m=1:M, s=1:S) for p=1:P]

    return EPP,  
           JuMP.objective_value(EPP), 
           JuMP.solve_time(EPP), 
           JuMP.relative_gap(EPP),
           JuMP.value.(E_toH2),
           JuMP.value.(E_fromH2),
           JuMP.value.(H_storage),
           JuMP.value.(η_toH2),
           JuMP.value.(η_fromH2),
           re,  # Return renewable energy data
           ED,  # Return energy demand
           e,   # Return energy demand matrix
           c_hBuy,  # Return H2 buy prices
           c_hSell; # Return H2 sell prices
end


# Refactor plotting function to use named parameters
function PlotResults(;alpha, X, E_toH2, E_fromH2, H_storage, E_charge, E_discharge, H_buy, H_sell, 
                    E_buy, E_sell, Periods, P=24, M=2, J=15, S=5, e=nothing, re=nothing,
                    c_hBuy=nothing, c_hSell=nothing)
    x = Periods

    # Graph 1: Scheduling (Machine States)
    function ablaufy1(m)
    y = [alpha[m,:,1] alpha[m,:,2] alpha[m,:,3] alpha[m,:,4] alpha[m,:,5]]*m
    end

    g1 = scatter(x, ablaufy1(1), label=["off" "turn on" "idle" "processing of job i" "turn off" "setup"], 
        markershape=[:diamond :diamond :utriangle :rect :diamond], colour=[:dodgerblue :red :orchid :limegreen :chocolate], markersize=10,
        ylabel="Machines",
        xlim=(0,P+1), ylim=(0.5,M+0.5), xticks=1:P, yticks=1:M, yflip=true, 
        left_margin=5Plots.mm, right_margin=5Plots.mm, top_margin=0Plots.mm, bottom_margin=0Plots.mm, size=(1000,500))

    for m=2:M
    scatter!(x,ablaufy1(m), label=:none, legend=:outerbottom, legend_font_pointsize=5,
            markershape=[:diamond :diamond :utriangle :rect :diamond], colour=[:dodgerblue :red :orchid :limegreen :chocolate], markersize=10)
    end

    for m=1:M, j=1:J, p=1:P
    annotate!((p, X[m,j,p]*m, text(string(j), :black, :10)));
    end
    
    # Common plot parameters
    plot_height = 500
    left_margin = 30Plots.mm

    # Graph 2a: Grid & Battery Energy Flows
    ED = zeros(length(x))
    if !isnothing(e)
        ED[x] = sum(alpha[m,x,s]*e[m,s] for m=1:M, s=1:S)
    end

    y_grid = [E_buy -E_sell E_charge E_discharge re -ED]
    g2a = groupedbar(x, y_grid, 
        label=["E_purchased" "E_sold" "E_charge" "E_discharge" "Erneuerbare" "Energiebedarf"],
        color=[:magenta3 :red2 :royalblue :lightblue :limegreen :darkred],
        ylabel="Energie [kWh]",
        title="Grid & Battery Energy Flows",
        xlim=(0,P), xticks=1:P,
        yticks=-200:50:200,
        left_margin=left_margin,
        right_margin=20Plots.mm,
        top_margin=5Plots.mm,
        bottom_margin=20Plots.mm,
        size=(1000,plot_height),
        legend=:outerbottom,
        legend_font_pointsize=8,
        framestyle=:box,
        grid=false,
        foreground_color_axis=nothing,
        foreground_color_border=nothing,
        tick_direction=:out)

    # Graph 2b: H2 System Energy Flows
    y_h2 = [-E_toH2 E_fromH2]  # Made E_toH2 negative to show below zero
    g2b = groupedbar(x, y_h2,
        label=["Elektrolyse" "Brennstoffzelle"],
        color=[:navy :skyblue],
        ylabel="Energie [kWh]",
        title="H₂ System Energy Flows",
        xlim=(0,P), xticks=1:P,
        yticks=-100:50:100,
        left_margin=left_margin,
        right_margin=20Plots.mm,
        top_margin=5Plots.mm,
        bottom_margin=20Plots.mm,
        size=(1000,plot_height),
        legend=:outerbottom,
        legend_font_pointsize=8,
        framestyle=:box,
        grid=false,
        foreground_color_axis=nothing,
        foreground_color_border=nothing,
        tick_direction=:out)

    # Graph 4: H2 Trading
    g4 = groupedbar(x, [H_buy -H_sell], 
        label=["H₂ Einkauf" "H₂ Verkauf"],
        color=[:green :red],
        ylabel="H₂ [kg]",
        title="H₂ Trading & Storage",
        xlim=(0,P), xticks=1:P,
        left_margin=left_margin,
        right_margin=20Plots.mm,
        top_margin=5Plots.mm,
        bottom_margin=20Plots.mm,
        size=(1000,plot_height),
        legend=:outerbottom,
        legend_font_pointsize=8,
        framestyle=:box,
        grid=false,
        foreground_color_axis=nothing,
        foreground_color_border=nothing,
        tick_direction=:out)

    # Graph 3: H2 Storage and Prices
    g3 = plot(x, H_storage, 
        label="H₂ Speicherstand",
        ylabel="H₂ [kg]", 
        color=:blue,
        title="H₂ Storage and Prices",
        xlim=(0,P), xticks=1:P,
        left_margin=left_margin,
        right_margin=20Plots.mm,
        top_margin=5Plots.mm,
        bottom_margin=20Plots.mm,
        size=(1000,plot_height),
        legend=:none,
        framestyle=:box,
        grid=false,
        foreground_color_axis=nothing,
        foreground_color_border=nothing,
        tick_direction=:out)
    
    # Add H2 prices on second y-axis
    plot!(twinx(), x, [c_hBuy c_hSell], 
        label=["H₂ Kaufpreis" "H₂ Verkaufpreis"],
        ylabel="Preis [€/kg]", 
        color=[:darkgreen :darkred], 
        line=:dash,
        legend=:none,
        framestyle=:box,
        grid=false,
        foreground_color_axis=nothing,
        foreground_color_border=nothing,
        tick_direction=:out)

    # Add combined legend for all three lines
    plot!(; legend=:outerbottom, 
          legend_font_pointsize=8,
          label=["H₂ Speicherstand" "H₂ Kaufpreis" "H₂ Verkaufpreis" "H₂ Einkauf" "H₂ Verkauf"])

    # Combine all plots
    plot(g1, g2a, g2b, g4, g3,
        layout=grid(5, 1), 
        size=(1200,2500),
        plot_title="Energy System Operation Schedule",
        plot_titlefontsize=14,
        plot_titlelocation=:center,
        top_margin=10Plots.mm,
        bottom_margin=20Plots.mm,
        subplot_spacing=0.05)
end

####################
# Main Code Section
####################

function main(directory_name::String="Test1", file_prefix::String="default")
    model, obj_val, solve_time, gap, E_toH2, E_fromH2, H_storage, η_toH2, η_fromH2, re, ED, e, c_hBuy, c_hSell = SolveEPP(360);
    
    H_sell = value.(model[:H_sell])
    H_buy = value.(model[:H_buy])
    
    results = Dict(
        "timestamp" => Dates.now(),
        "objective_value" => obj_val,
        "solve_time" => solve_time,
        "gap" => gap,
        "avg_h2_production_efficiency" => mean(η_toH2),
        "avg_h2_usage_efficiency" => mean(η_fromH2),
        "total_E_toH2" => sum(E_toH2),
        "total_E_fromH2" => sum(E_fromH2),
        "total_H_storage" => sum(H_storage),
        "total_H_sell" => sum(H_sell),
        "total_H_buy" => sum(H_buy)
    )
    
    detailed_results = DataFrame(
        period = 1:24,
        renewable_energy = re,
        energy_demand = ED,
        E_toH2 = E_toH2,
        E_fromH2 = E_fromH2,
        H_storage = H_storage,
        H_sell = H_sell,
        H_buy = H_buy,
        η_toH2 = η_toH2,
        η_fromH2 = η_fromH2
    )

    results_dir = "results/$(directory_name)"
    mkpath(results_dir)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    summary_file = "$(results_dir)/$(file_prefix)_summary_$(timestamp).txt"
    detailed_file = "$(results_dir)/$(file_prefix)_detailed_$(timestamp).csv"
    plot_file = "$(results_dir)/$(file_prefix)_plot_$(timestamp).png"
    
    open(summary_file, "w") do io
        println(io, "#########################################################################################")
        println(io, "OPTIMIZATION RESULTS:")
        println(io, "Timestamp: ", results["timestamp"])
        println(io, "Objective value: ", results["objective_value"])
        println(io, "Solve time: ", results["solve_time"], " seconds")
        println(io, "Relative gap: ", results["gap"])
        println(io, "\nEFFICIENCY METRICS:")
        println(io, "H2 production efficiencies by period: ", η_toH2)
        println(io, "H2 usage efficiencies by period: ", η_fromH2)
        println(io, "\nPERIOD-BY-PERIOD VALUES:")
        println(io, "\n" * "#"^150)  
        @printf(io, "%-6s | %12s | %12s | %12s | %12s | %12s | %12s | %12s\n", 
               "Period", "E_toH2", "E_fromH2", "H_storage", "H_sell", "H_buy", "E_charge", "E_discharge")
        println(io, "-"^150)  
        for p in 1:24
            @printf(io, "%-6d | %12.2f | %12.2f | %12.2f | %12.2f | %12.2f | %12.2f | %12.2f\n", 
                   p, E_toH2[p], E_fromH2[p], H_storage[p], H_sell[p], H_buy[p], value(model[:E_charge][p]), value(model[:E_discharge][p]))
        end
        println(io, "#"^150)  
        println(io, "\nSUMMARY STATISTICS:")
        println(io, "Average H2 production efficiency: ", results["avg_h2_production_efficiency"])
        println(io, "Average H2 usage efficiency: ", results["avg_h2_usage_efficiency"])
        println(io, "Total E_toH2: ", results["total_E_toH2"])
        println(io, "Total E_fromH2: ", results["total_E_fromH2"])
        println(io, "Total H_storage: ", results["total_H_storage"])
        println(io, "Total H_sell: ", results["total_H_sell"])
        println(io, "Total H_buy: ", results["total_H_buy"])
        println(io, "#########################################################################################")
    end
    
    CSV.write(detailed_file, detailed_results)
    
    println("starting test run ...")
    println("#########################################################################################")
    println("zusammenfassung:")
    println("Objective value: ", obj_val, " found after ", solve_time, " seconds. Relative gap is: ", gap)
    println("")
    println("Creating visualization...")
    p = PlotResults(
        alpha=value.(model[:alpha]),
        X=value.(model[:X]),
        E_toH2=E_toH2,
        E_fromH2=E_fromH2,
        H_storage=H_storage,
        E_charge=value.(model[:E_charge]),
        E_discharge=value.(model[:E_discharge]),
        H_buy=H_buy,
        H_sell=H_sell,
        E_buy=value.(model[:E_buy]),
        E_sell=value.(model[:E_sell]),
        Periods=1:24,
        e=e,
        re=re,
        S=5,
        c_hBuy=c_hBuy,
        c_hSell=c_hSell
    );
    savefig(p, plot_file)
    println("Results saved in $(results_dir)/ directory with timestamp $(timestamp)")
    println()
end

# wenn man die datai ausführt soll main ausgeführt werden
if abspath(PROGRAM_FILE) == @__FILE__
    main("Test7", "neueChargeDischargeKostn")
end

#main("MyDirectory", "MyTest"); # Uncomment and modify to run with custom names