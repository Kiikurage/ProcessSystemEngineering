//
//  main.cpp
//  work
//
//  Created by KikuraYuichiro on 2014/10/22.
//  Copyright (c) 2014年 KikuraYuichiro. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <math.h>

#define PRODUCT_FLOW_RATE__             63.0    //  Gp[ton/day]
#define TOTAL_CONVERSION_RATE           0.6     //  ζ
#define RESIDENCE_TIME_IN_ONE_REACTOR__ 5       //  τ [hour]
#define MIXTURE_DENSITY                 850.0   //  ρ [kg/m^3]
#define ASPECT_RATE                     1.3     //  α
#define POLYMERIZATION_HEAT__           72.8    //  ⊿H[kJ/mol]
#define MIXTURE_TEMP                    50.0    //  T1[℃]
#define TOLUENE_THERMAL_CONDUCTIVITY    0.128   //  λ [W/K*m]
#define TOLUENE_SPECIFIC_HEAT__         1.68    //  Cp[J/g*K]
#define TOLUENE_MOLECULAR_WEIGHT__      54.0    //  M [g/mol]
#define ML                              50      //  Monomer Length

//unit-converted parameters
#define PRODUCT_FLOW_RATE               (PRODUCT_FLOW_RATE__*1.0e3/24/60/60)    //  Gp[kg/s]
#define RESIDENCE_TIME_IN_ONE_REACTOR   (RESIDENCE_TIME_IN_ONE_REACTOR__*60*60) //  τ0[s]
#define POLYMERIZATION_HEAT             (POLYMERIZATION_HEAT__*1.0e3)           //  ⊿H[J/mol]
#define TOLUENE_SPECIFIC_HEAT           (TOLUENE_SPECIFIC_HEAT__*1.0e3)         //  Cp[J/kg*K]
#define TOLUENE_MOLECULAR_WEIGHT        (TOLUENE_MOLECULAR_WEIGHT__*1.0e-3)     //  M [kg/mol]

#define printLine printf("\n----------------------------------------------------------------\n\n")

float getPowerConsumption(int reactorNumber,
                          double initialConcentration,
                          double coolantTemp,
                          bool flagLog) {
    
    if (flagLog) {
        printLine;
        printf("[Input Parameter]\n\n");
        printf("  number of reactor                : N  = %6d\n",           reactorNumber);
        printf("  initial concentration of monomer : γ0 = %10.3lf\n",        initialConcentration);
        printf("  temperature of coolant           : T2 = %10.3lf[°C]\n",    coolantTemp);
    }

    //----------------------------------------------------------------
    //  calculate common reactor information.

    double z = 1.0/(1.0 - TOTAL_CONVERSION_RATE);
    double reactionRateConst = 1.0 / RESIDENCE_TIME_IN_ONE_REACTOR * (z - 1);
    double residenceTime = (pow(z, 1.0/reactorNumber) - 1) / (z - 1) * RESIDENCE_TIME_IN_ONE_REACTOR;
    double totalResidenceTime = residenceTime * reactorNumber;

    double feed = PRODUCT_FLOW_RATE / (MIXTURE_DENSITY * initialConcentration * TOTAL_CONVERSION_RATE);
    double volume = feed * residenceTime;
    
    double diameter = cbrt(4.0*volume/ASPECT_RATE/M_PI);
    double height = diameter * ASPECT_RATE;
    double sideArea = M_PI*diameter*height;
    
    //----------------------------------------------------------------
    //  output common reactor information.
    if (flagLog) {
        printLine;
        printf("[Constances]\n");
        printf("  product flow rate                : Gp = %10.3lf[ton/day]\n",   PRODUCT_FLOW_RATE__);
        printf("  total conversion rate            : ζ  = %10.3lf\n",            TOTAL_CONVERSION_RATE);
        printf("  residence time of one reactor    : τ  = %10.3lf[s]\n",         residenceTime);
        printf("  mixture density                  : ρ  = %10.3lf[kg/m^3]\n",    MIXTURE_DENSITY);
        printf("  reactor aspect rate              : α  = %10.3lf\n",            ASPECT_RATE);
        printf("  polymerization heat              : ⊿H = %10.3lf[kJ/mol]\n",    POLYMERIZATION_HEAT__);
        printf("  mixture temperature              : T1 = %10.3lf[°C]\n",        MIXTURE_TEMP);

        printf("  thermal conductivity of toluene  : λ  = %10.3lf[W/K*m]\n",     TOLUENE_THERMAL_CONDUCTIVITY);
        printf("  specific heat of toluene         : Cp = %10.3lf[J/K*g]\n",     TOLUENE_SPECIFIC_HEAT__);
        printf("  molecular weight of toluene      : M  = %10.3lf[g/mol]\n",     TOLUENE_MOLECULAR_WEIGHT__);
        printf("  monomer length                   : ML = %6d\n",                ML);
        printf("  reaction rate const              : K  = %10.3lf[1/s]\n",       reactionRateConst);
        
        printLine;
        printf("[Reactors Overview]\n");
        printf("  feed rate                        : F  = %10.3f[m^3/s]\n",      feed);
        printf("  volume                           : V  = %10.3f[m^3]\n",        volume);
        printf("  diameter                         : D  = %10.3lf[m]\n",         diameter);
        printf("  height                           : H  = %10.3lf[m]\n",         height);
        printf("  side area                        : S  = %10.3lf[m^2]\n",       sideArea);
    }
    
    //----------------------------------------------------------------
    //  calculate for each reactant value
    
    int i;
    double
        concentRation,
        conversionRate,
        productHeat,
        heatTransCof,
        viscosity,
        Nu,
        Pr,
        Re,
        Np,
        revolutionNumber,
        powerConsumption,
        totalPowerConsumption = 0;
    
    for (i = 1; i <= reactorNumber; i++) {
        concentRation = initialConcentration * pow(z, -1.0*i/reactorNumber);
        
        productHeat = POLYMERIZATION_HEAT * MIXTURE_DENSITY * concentRation / TOLUENE_MOLECULAR_WEIGHT * reactionRateConst * volume;
        
        heatTransCof = productHeat / sideArea / (MIXTURE_TEMP - coolantTemp);
        
        conversionRate = 1.0 - pow(1.0 - TOTAL_CONVERSION_RATE, 1.0*i/reactorNumber);

        viscosity = pow(ML, 1.7) * pow(conversionRate, 2.5) * exp(21 * initialConcentration) * 1.0e-3;
        
        Nu = heatTransCof * diameter / TOLUENE_THERMAL_CONDUCTIVITY;
        Pr = viscosity / TOLUENE_THERMAL_CONDUCTIVITY * TOLUENE_SPECIFIC_HEAT;
        Re = pow(2.0*Nu*pow(Pr, -1.0/3), 3.0/2);
        Np = 14.6 * pow(Re, -0.28);
        revolutionNumber = Re * ((viscosity / MIXTURE_DENSITY) / pow(diameter/2, 2));
        
        powerConsumption = Np * MIXTURE_DENSITY * pow(revolutionNumber, 3)*pow(diameter/2, 5);
        
        totalPowerConsumption += powerConsumption;

        //----------------------------------------------------------------
        //  output for each reactant value
        
        if (flagLog) {
            printLine;
            printf("[Reactor %d]\n",                                                i);
            printf("  concentration of Butadien of out : C_%-2d = %10.3lf\n",       i,  concentRation);
            printf("  product heat of polymerization   : Qp   = %10.3lf[W]\n",          productHeat);
            printf("  viscosity                        : μ_%-2d = %10.3lf[Pa*s]\n", i,  viscosity);
            printf("  Nusselt number                   : Nu   = %10.3lf\n",             Nu);
            printf("  Prandtl number                   : Pr   = %10.3lf\n",             Pr);
            printf("  Reynolds number                  : Re   = %10.3lf\n",             Re);
            printf("  Power number                     : Np   = %10.3lf\n",             Np);
            printf("  revolution number                : n_%-2d = %10.3lf[rps]\n",  i,  revolutionNumber);
            printf("  power consumption                : P_%-2d = %10.3lf[W]\n",    i,  powerConsumption);
        }
    }
    
    //----------------------------------------------------------------
    //  output for each reactant value
    
    if (flagLog) {
        printLine;
        printf("[Total]\n");
        printf("  total power consumption          : P  = %10.3f[W]\n", totalPowerConsumption);
    }

    if (flagLog) {
        printLine;
    }

    return totalResidenceTime;
}

void showUsage(const char *argv[]) {
    fprintf(stderr, "Usage : %s N y0 T2 [--log]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    N : number of reactor\n");
    fprintf(stderr, "    γ0: initial concentration of monomer\n");
    fprintf(stderr, "    T2: temperature of coolant [°C]\n");
    fprintf(stderr, " --log: when true, output full logs.\n");
    fprintf(stderr, "\n");
}

int main(int argc, const char *argv[]) {
    int reactorNumber;
    double initialConcentration, coolantTemp;
    bool flagLog = false;
    
    flagLog = (argc == 5 && (strcmp(argv[4], "--log") == 0));
    
    if (!((argc == 4) || (flagLog && argc == 5))) {
        showUsage(argv);
        exit(-1);
    }

    reactorNumber = atoi(argv[1]);
    initialConcentration = atof(argv[2]);
    coolantTemp = atof(argv[3]);
    
    printf("%lf\n", getPowerConsumption(reactorNumber, initialConcentration, coolantTemp, flagLog));

    return 0;
}