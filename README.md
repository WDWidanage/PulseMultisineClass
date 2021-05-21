# Pulse Multisine Class - To estimate non-linear equivalent circuit models of lithium-ion batteries


## Intro
This MATALB class allows a user to define a pulse-multisine signal (as an alternative to the commonly used HPPC - Hybrid Pulse Power Characterisation test) to develop a non-linear equivalent circuit model of lithium-ion batteries as a function of state-of-charge and temperature.

## Step 1 - Generate a pulse-multisine
The procedure to design a pulse-multisine signal is explained in [Paper 1]( https://www.sciencedirect.com/science/article/pii/S0378775316305511?casa_token=j3psTFsv61QAAAAA:R2WsiBDDyz6AVdILwqrGIv6TMKm21G8d2I-9FcbzAJIKShRsLEp6GjY1GGoNyTtQjsdxLVdmoQ)

Clone or download the repository and run “Generate_PulseMultisine.m”. The steps and further explanations are provided in the doc string of the script.

## Step 2 - Estimate the NL-ECM model
The procedure to estimate a non-linear equivalent circuit model is explained in [Paper 2]( https://www.sciencedirect.com/science/article/pii/S037877531630550X?casa_token=c9NRYDoaX00AAAAA:2hul3zvGP6A9JvKiHDUKTxpE_qjb5zY3guZ_UCpUYPcCHcYvDPN9_QiCHQB5wXYDhyfbNE8rew)

Clone or download the repository and run “estimateNLECM.m”. The steps and further explanations are provided in the doc string of the script.

## Step 3 - Enjoy and have a cuppa!
