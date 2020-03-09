### Domo de Saturacao ###

# para obtencao do domo sao necessarios ℘ V1 V2 vr1 vr2 Tr_sat Pr_sat

#V1(vr1) = 3vr1 - 1

#V2(vr2) = 3vr2 - 1

#℘(V1,V2) = V1 / V2

#V1(℘) = ((℘^2)*(-2℘ + (℘ + 1)*log(℘) + 2)) / ((℘^2) - 2*℘*log(℘) -1)

#V2(℘) = (℘*(-2℘ + (℘ + 1)*log(℘) + 2)) / ((℘^2) - 2*℘*log(℘) -1)

V2(℘) = (2*log(℘)/(℘-1) - 1/℘ - 1) / (2 - (℘+1)*log(℘)/(℘-1))

V1(℘) = V2(℘)*℘

vr1(V1) = (V1 + 1) / 3

vr2(V2) = (V2 + 1) / 3

Tr_sat(V1,V2) = (27*V1*V2*(V1 + V2 + 2)) / (8*((V1 + 1)^2)*((V2 + 1)^2))

Pr_sat(V1,V2) = (27*(-V1 - V2 + (V1 + 1)*(V2 + 1) - 2)) / (((V1 + 1)^2)*((V2 + 1)^2))

℘list = range(0, stop =0.95, length = 20)

V2list = []

V1list = []

vr1list = []

vr2list = []

