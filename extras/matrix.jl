c_hBuy = [5.74519219, 7.86604509, 6.73319315, 7.05671601, 6.74685032, 5.22883979,
6.83608179, 8.46972794, 5.03387753, 4.53276037, 7.03822001, 5.89296313,
3.96507793, 6.05017275, 4.2303096, 4.80617114, 5.63683564, 6.45251204,
3.70139256, 6.59991497, 4.95218885, 7.20329569, 4.75833894, 5.81176465] * 1;


c_hSell = [3.44711531, 4.71962705, 4.03991589, 4.2340296, 4.04811019, 3.13730387,
4.10164907, 5.08183676, 3.02032652, 2.71965622, 4.222932, 3.53577788,
2.37904676, 3.63010365, 2.53818576, 2.88370269, 3.38210138, 3.87150722,
2.22083554, 3.95994898, 2.97131331, 4.32197742, 2.85500336, 3.48705879] * 1


x = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
c_sell = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5] *0.02;
println(c_sell)

length(x)

365 * 24

300 * 1/33.3333
# euro / kwh * 33.3333333 kwh / kg
0.15 * 33.3333333

1 * 1/0.6 


rounded_c_hBuy = round.(c_hBuy)
rounded_c_hSell = round.(c_hSell)


println("c_hBuy = ", rounded_c_hBuy)
println("c_hSell = ", rounded_c_hSell)


leistung_punkte              = [0,   10,    20,    30,    40,    50,   60,    70,   80,   90,   100]   # kW (Leistung)
wirkungsgrad_punkte_toH2     = [0, 0.95,  0.90,  0.88,  0.85,  0.82, 0.80,  0.78, 0.75, 0.72,  0.70]   # toH2
wirkungsgrad_punkte_fromH2   = [0,  0.6,  0.58,  0.55,  0.52,  0.50, 0.48,  0.45, 0.43,  0.4,   0.3]   # fromH2

1/33.333333

1/3

100 * 0.3 * 1/33.3333
100 * 0.8 * 1/33.3333


E_toH2 * 

@constraint(EPP, [p=1:1], H_storage[p] == H_storage_initial + EC_menge[p] - FC_menge[p] + H_buy[p] - H_sell[p])
# other periods
@constraint(EPP, [p=2:P], H_storage[p] == H_storage[p-1] + EC_menge[p] - FC_menge[p] + H_buy[p] - H_sell[p])
# last period
@constraint(EPP, [p=P:P], H_storage_end == H_storage[p-1] + EC_menge[p] - FC_menge[p] + H_buy[p] - H_sell[p])

@constraint(EPP, [p=1:P], EC_menge[p] == E_toH2[p] * η_toH2[p] * H_heizwert)  # kWh * % * kg/kwh = kg
@constraint(EPP, [p=1:P], FC_menge[p] == E_fromH2[p] * 1/η_fromH2[p] * H_heizwert) 

# 0 kw -> 0.0001%  
1/0
function safe_inverse(x)
    return x == 0 ? 0 : 1/x  
end
wirkungsgrad_punkte_toH2_inv = [safe_inverse(x) for x in wirkungsgrad_punkte_toH2]
wirkungsgrad_punkte_fromH2_inv = [safe_inverse(x) for x in wirkungsgrad_punkte_fromH2]