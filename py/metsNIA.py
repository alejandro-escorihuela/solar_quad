#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# 05-06-2024
# alex
# metsNIA.py

import numpy as np
from extensio import *

## Escissió general ABA
SVABA = ([0.5, 0.5], [1.0]) ###### ja està estes

## Escissió RKN ABA
NA7 = [0.11461549414048543, 0.5342629313313436, 0.17711460144207827, 0.3777694345548009, 0.29421383706768817, -0.19724296168040745]
NA14 = [0.0378593198406116, 0.102635633102435, -0.0258678882665587, 0.314241403071447, -0.130144459517415, 0.106417700369543, -0.00879424312851058, 0.09171915262446165, 0.183983170005006, -0.05653436583288827, 0.004914688774712854, 0.143761127168358, 0.328567693746804]
NA19 = list(np.longfloat(["0.0505805", "0.149999", "-0.0551795510771615573511026950361", "0.423755898835337951482264998051",
                          "-0.213495353584659048059672194633", "-0.0680769774574032619111630736274", "0.227917056974013435948887201671", "-0.235373619381058906524740047732",
                          "0.387413869179878047816794031058",
                          "0.129478606560536730662493794395", "0.222257260092671143423043559581", "-0.0577514893325147204757023246320", "-0.0578312262103924910221345032763",
                          "0.103087297437175356747933252265", "-0.140819612554090768205554103887", "0.0234462603492826276699713718626", "0.134854517356684096617882205068",
                          "0.0287973821073779306345172160211"]))

## NIA
NIA2_42 = ([(3-np.sqrt(3))/6, 1.0 - (3-np.sqrt(3))/3, (3-np.sqrt(3))/6], [0.5, 0.5]) ###### ja està estes
NIA6_122 = [0.033765242898423986093849222753002695, 0.135630063868443757075450979737044631, 0.211295100191533802515448936669596706, 0.085662246189585172520148071086366447, 0.180380786524069303784916756918858056]
NIA7_864 = [0.0711334264982231177779387300061549964174, 0.241153427956640098736487795326289649618, 0.521411761772814789212136078067994229991, 0.183083687472197221961703757166430291072, 0.310782859898574869507522291054262796375, -0.0265646185119588006972121379164987592663]
NIA8_1064 = [0.0380944974224122, 0.1452987161169130, 0.2076276957255412, 0.4359097036515262, 0.0958588808370752, 0.2044461531429988, 0.2170703479789911]

## NIA nous
NIA884 = [0.06669891497364668, 0.2359095515194927, -0.10626513574360669, -0.16666163973367262, -0.017886377696568687, 0.16644076324010687, 0.4278213876724725, -0.1253625358648917, 0.12279013970528008]
NIA1084 = [0.06458266253564669, 0.2628461314694317, -0.23936775641840224, -0.01253197543793027, 0.1689208508190114, 0.16258030603206056, -0.16307005840850436, -0.3337680728967171, 0.26274501629762026, 0.38278068984047986]
NIA1086A = [0.2618013676833966, -0.15370001231853378, -0.0493509922784518, -0.020840998612185847, 0.3058195065158628, 0.23433789395941815, -0.013569437747433345, 0.28762803979298224, -0.2991984339495876, 0.23425735986406734, 0.2924869991673998, -0.1062469892771738]
NIA1086B = [0.12709234820095489, -0.06025822519628227, 0.13554447235067524, -0.15759421307537005, 0.3232455982705558, 0.17437198974126544, 0.2886475755274505, -0.28581393064631744, -0.024020497448570086, 0.2520639565498974, 0.3168491397542202, -0.26043799397390255]
## RNIA nous
RNIA1086 = [0.03742418824257127, 0.2656106460487173, -0.11898009697774081, 0.25934295597727025, -0.05048710950481768, -0.17579709169189453, 0.09462046562180979, -0.08311300330604939, 0.1708689632198862, -0.09852232690976058, 0.3704703757050709]

## Nucli NIA
NIA3_x64 = [0.5600879810924619, 1.5171479707207228]

## Starter NIA
TNIA8_A = [0.05490386371597787, -0.28283514951371963, 1.3912359747629677, -0.31031471187617554, -0.17411209690059248, 0.32112211981154207999, 0.2924509441030227, -0.019451859711026816, 0.036315562429046384, -0.29203995279233447, 1.035567080563363, -0.05284177459207079798]
TNIA8_B = [0.015855169600235627, 0.5760860447577043, -1.0454628537953747, 1.4535774759561593, 1.173891649242811, -1.1739474857615355269, 0.21998028005693096, 0.8031998140452785, -0.0064852633141201975, -28.088878551868973, 0.0001688373612304016, 28.072014883719653335]
TNIA10 = [0.13854271338364255, 1.7638286294297816, -0.02226826621201071, -0.4153821320829431, 0.23971791586543062, -1.0286391118497091, 0.32420025146580814005, 0.34953719460449373, 0.9721369343176093, -1.25500602179785, -0.2919356517110871, 0.49539466687406475, 0.6819181447433431, 0.04795473296942621993]

## Processat NIA
PNIA14 = [-0.2822177053805228, -0.18331921985438093, 0.7832862612309347, 0.07693317031111557, 0.07019285265291131, -0.3832728226411798, -0.21418604065030244, -0.36180460311212004, 0.4818528837549031, 0.1483521570725569, -0.4739645520570135, 0.5322008825362824, -0.10539686980026997, 0.3735875862057152, -0.8373460255853887, 0.501599467256724, -0.47909755935255793, 0.8501031588160132, -0.7804811294241945, 0.2829781080207743, -1.2169214003761924e-05, 0.029196376454051577, -0.06083582987887891, 0.04740565477100147, 0.0030133129557245613, -0.06339107910207395, 0.023703524098908944, -0.011451160188582951, -0.26679476605855695, -0.04844467342082413, 0.14022836433710442, 0.0924088212990693, 0.0023038305565269185, -0.02742838135009622, -0.15813236603903924, -0.007110515476561043, 0.0761744117645311, 0.007663603614189092, -0.08936298894581002, 0.3108660298233198]


# Mètodes a emprar
cofABA, proABA = [], []
cofABA.append(estendreA(NIA7_864))
proABA.append([])
cofABA.append(estendreA(NIA8_1064))
proABA.append([])    
# cofABA.append(estendreA(NIA8_1064))
# proABA.append(PNIA14)
# cofABA.append(estendreA(NIA884))
# proABA.append([])
# cofABA.append(estendreA(NIA1084))
# proABA.append([])
cofABA.append(estendreA(NIA1086A))
proABA.append([])
cofABA.append(estendreA(RNIA1086))
proABA.append([])
cofABA.append(estendreA(NIA3_x64))
proABA.append(TNIA8_A)
# cofABA.append(estendreA(NIA3_x64))
# proABA.append(TNIA8_B)
# cofABA.append(estendreA(NIA3_x64))
# proABA.append(TNIA10)     
# cofABA.append(estendreA(NA7))
# proABA.append([])
# cofABA.append(estendreA(NA19))
# proABA.append([])     

noms = ["NIA7_864", "NIA8_1064", "NIA13_1086", "RNIA12_1086", "TNIA3_864"]
