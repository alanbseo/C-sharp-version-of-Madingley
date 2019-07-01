using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

namespace Madingley
{
    /// <summary>
    /// Removes autotroph matter appropriated by humans from a grid cell's autotroph stocks
    /// </summary>
    /// <remarks>Assumes that autotroph matter is appropriated evenly from different stocks in proportion to their biomass</remarks>
    public class HumanAutotrophMatterAppropriation
    {
        //Global mean HANPP from harvesting values from Haberl et al. (2007), PNAS 
        // Table 2. Breakdown of global HANPP harvest (excluding human-induced fires) in the year 2000 to land-use classes
        // gC/m2/yr
        double HANPPh_crop = 296;
        double HANPPh_wilderness = 0.0;
        double HANPPh_grazing = 41;
        double HANPPh_forestry = 48;
        double HANPPh_urban = 63;

        UtilityFunctions _Utilities;

        /// <summary>
        /// Constructor for human appropriation of autotroph matter
        /// </summary>
        public HumanAutotrophMatterAppropriation()
        {

            _Utilities = new UtilityFunctions();
        }

        /// <summary>
        /// Remove human appropriated matter from the grid cell autotroph stocks
        /// </summary>
        /// <param name="cellEnvironment">The environment in the current grid cell</param>
        /// <param name="humanNPPScenario">The type of NPP extraction to apply</param>
        /// <param name="gridCellStocks">The stocks in the current grid cell</param>
        /// <param name="actingStock">The position of the acting stock in the jagged array of grid cell stocks</param>
        /// <param name="currentTimestep">The current model time step</param>
        /// <param name="burninSteps">The number of steps to run before impact is simulated</param>
        /// <param name="impactSteps">The number of time steps to apply the impact for (for 'temporary' scenarios)</param>
        /// <param name="impactCell">Whether this cell should have human impacts applied</param>
        /// <remarks>Scenario types are: 'no' = no removal; 'hanpp' = appropriated NPP estimate from input map; constant = constant appropriation after burn-in; 
        /// temporary = constant after burn-in until specified time; value = proportion of plant biomass appropriated</remarks>
        public double RemoveHumanAppropriatedMatter(double wetMatterNPP, SortedList<string, double[]> cellEnvironment,
            Tuple<string, double, double> humanNPPScenario, GridCellStockHandler
            gridCellStocks, int[] actingStock, uint currentTimestep, int scenarioYear, uint burninSteps,
            uint impactSteps, uint recoverySteps, uint instantStep, uint numInstantStep, Boolean impactCell,
            string globalModelTimestepUnits, FunctionalGroupDefinitions madingleyStockDefinitions, double fracEvergreen)
        {

            double RemovalRate = 0.0;

            if (impactCell)
            {
                // Factor to convert NPP from units per m2 to units per km2
                double m2Tokm2Conversion = 1000000.0;


                if (humanNPPScenario.Item1 == "hanpp") //static hanpp 
                {

                    if (currentTimestep > burninSteps)
                    {
                        // Loop over stocks in the grid cell and calculate the total biomass of all stocks
                        double TotalAutotrophBiomass = 0.0;
                        foreach (var stockFunctionalGroup in gridCellStocks)
                        {
                            for (int i = 0; i < stockFunctionalGroup.Count; i++)
                            {
                                TotalAutotrophBiomass += stockFunctionalGroup[i].TotalBiomass;
                            }
                        }

                        // Get the total amount of NPP appropriated by humans from this cell
                        double HANPP = cellEnvironment["HANPP"][0];

                        // If HANPP value is missing, then assume zero
                        if (HANPP == cellEnvironment["Missing Value"][0]) HANPP = 0.0;

                        HANPP *= cellEnvironment["Seasonality"][currentTimestep % 12];

                        // Allocate HANPP for this stock according to the proportion of total autotroph biomass that the stock represents
                        if (TotalAutotrophBiomass == 0.0)
                        {
                            HANPP = 0.0;
                        }
                        else
                        {
                            HANPP *= (gridCellStocks[actingStock].TotalBiomass / TotalAutotrophBiomass);
                        }


                        // Convert gC/m2/month to gC/km2/month
                        HANPP *= m2Tokm2Conversion;

                        // Multiply by cell area (in km2) to get g/cell/month
                        HANPP *= cellEnvironment["Cell Area"][0];


                        // Convert from gC to g dry matter
                        double DryMatterAppropriated = HANPP * 2;

                        // Convert from g dry matter to g wet matter
                        double WetMatterAppropriated = DryMatterAppropriated * 2;


                        //Calculate the rate of HANPP offtake
                        if (wetMatterNPP.CompareTo(0.0) == 0)
                        {
                            RemovalRate = 0.0;
                        }
                        else
                        {
                            RemovalRate = Math.Min(1.0, WetMatterAppropriated / wetMatterNPP);
                        }
                        // Remove human appropriated autotroph biomass from total autotroph biomass
                        //gridCellStocks[actingStock].TotalBiomass -= WetMatterAppropriated;

                        //if (gridCellStocks[actingStock].TotalBiomass < 0.0) gridCellStocks[actingStock].TotalBiomass = 0.0;
                    }
                }
                else if(humanNPPScenario.Item1 == "ssp") // annual hanpp input
                {
                    //The scenario year calculation removes the need for this if check and allows the burnin period to have HANPP applied
                    //if (currentTimestep > burninSteps)
                    {
                        // Get the total amount of NPP appropriated by humans from this cell
                        double HANPPh = cellEnvironment["HANPPharvest"][scenarioYear];
                        double HANPPlc = cellEnvironment["HANPPlc"][scenarioYear];

                        // If HANPP value is missing, then assume zero
                        if (HANPPh == cellEnvironment["Missing Value"][0]) HANPPh = 0.0;
                        if (HANPPlc == cellEnvironment["Missing Value"][0]) HANPPlc = 0.0;

                        HANPPh *= cellEnvironment["Seasonality"][currentTimestep % 12];
                        HANPPlc *= cellEnvironment["Seasonality"][currentTimestep % 12];
                        
                        //Allocate HANPP to this stock depending on its definition
                        
                        //Allocate HANPP between deciduous and evergreen according to FracEvergreen
                        if (madingleyStockDefinitions.GetTraitNames("leaf strategy", actingStock[0]).Equals("deciduous"))
                        {
                            HANPPlc *= (1.0 -fracEvergreen);
                            HANPPh *= (1.0 - fracEvergreen);
                        }
                        else
                        {
                            HANPPlc *= (fracEvergreen);
                            HANPPh *= (fracEvergreen);
                        }
                        //Allocate between impacted and natural using fractional cell area
                        if (madingleyStockDefinitions.GetTraitNames("impact state",actingStock[0]).Equals("primary"))
                        {
                            //No Hanpp from land cover change in natural lands
                            HANPPlc *= 0.0;
                            HANPPh *= 0.0;
                        } else if (madingleyStockDefinitions.GetTraitNames("impact state", actingStock[0]).Equals("secondary"))
                        {
                            //No Hanpp from land cover change in natural lands
                            HANPPlc *= 0.0;
                            HANPPh *= 1.0 - FracImpactedHANPPh(cellEnvironment["Fsecondary"][scenarioYear], cellEnvironment["Fcropland"][scenarioYear],
                                cellEnvironment["Furban"][scenarioYear], cellEnvironment["Fgrazing"][scenarioYear]);
                        }
                        else
                        {
                            //All HANPPlc comes from impacted lands
                            HANPPh *= (FracImpactedHANPPh(cellEnvironment["Fsecondary"][scenarioYear], cellEnvironment["Fcropland"][scenarioYear],
                                cellEnvironment["Furban"][scenarioYear], cellEnvironment["Fgrazing"][scenarioYear]));
                        }

                        //Combine harvest and land change terms
                        double HANPP = HANPPh + HANPPlc;

                        // Convert gC/m2/month to gC/km2/month
                        HANPP *= m2Tokm2Conversion;

                        // Multiply by cell area (in km2) to get g/cell/month
                        HANPP *= gridCellStocks[actingStock].FractionalArea * cellEnvironment["Cell Area"][0];


                        // Convert from gC to g dry matter
                        double DryMatterAppropriated = HANPP * 2;

                        // Convert from g dry matter to g wet matter
                        double WetMatterAppropriated = DryMatterAppropriated * 2;


                        //Calculate the rate of HANPP offtake (%)
                        if (wetMatterNPP.CompareTo(0.0) == 0)
                        {
                            RemovalRate = 0.0;
                        }
                        else
                        {
                            RemovalRate = Math.Min(1.0, WetMatterAppropriated / wetMatterNPP);
                        }

                    }


                }
                else if (humanNPPScenario.Item1 == "no")
                {
                    // Do not remove any autotroph biomass
                }
                else if (humanNPPScenario.Item1 == "constant")
                {
                    // If the burn-in period has been completed, then remove the specified constant
                    // fraction from the acting autotroph stock
                    if (currentTimestep > burninSteps)
                    {
                        //gridCellStocks[actingStock].TotalBiomass -= (gridCellStocks[actingStock].TotalBiomass *
                        //    humanNPPScenario.Item2);
                        RemovalRate = humanNPPScenario.Item2;
                    }
                }
                else if (humanNPPScenario.Item1 == "temporary")
                {
                    // If the spin-up period has been completed and the period of impact has not elapsed,
                    // then remove the specified constant fraction from the acting autotroph stock
                    if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps)))
                    {
                        //gridCellStocks[actingStock].TotalBiomass -= (gridCellStocks[actingStock].TotalBiomass *
                        //    humanNPPScenario.Item2);
                        RemovalRate = humanNPPScenario.Item2;
                    }

                }
                else if (humanNPPScenario.Item1 == "escalating")
                {
                    // If the spin-up period has been completed, then remove a proportion of plant matter
                    // according to the number of time-steps that have elapsed since the spin-up ended
                    if (currentTimestep > burninSteps)
                    {
                        //gridCellStocks[actingStock].TotalBiomass -= gridCellStocks[actingStock].TotalBiomass *
                        //    (Math.Min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2)));

                        RemovalRate = (Math.Min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2)));

                    }
                }
                else if (humanNPPScenario.Item1 == "temp-escalating")
                {
                    // If the spin-up period has been completed and the period of impact has not elapsed, 
                    // then remove a proportion of plant matter
                    // according to the number of time-steps that have elapsed since the spin-up ended
                    if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps)))
                    {
                        //gridCellStocks[actingStock].TotalBiomass -= gridCellStocks[actingStock].TotalBiomass *
                        //    (Math.Min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2)));

                        RemovalRate = (Math.Min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2)));
                    }

                }
                else if (humanNPPScenario.Item1 == "temp-escalating-const-rate")
                {
                    // If the spin-up period has been completed and the period of impact (specified by the third scenario element
                    // has not elapsed, 
                    // then remove a proportion of plant matter
                    // according to the number of time-steps that have elapsed since the spin-up ended

                    int ConstImpactSteps = Convert.ToInt32(humanNPPScenario.Item3 * _Utilities.ConvertTimeUnits("year", globalModelTimestepUnits));

                    if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + ConstImpactSteps)))
                    {
                        //gridCellStocks[actingStock].TotalBiomass -= gridCellStocks[actingStock].TotalBiomass *
                        //    (Math.Min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2)));

                        RemovalRate = (Math.Min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2)));
                    }
                }
                else if (humanNPPScenario.Item1 == "temp-escalating-const-rate-duration")
                {
                    // If the spin-up period has been completed and the period of impact (specified by the third scenario element
                    // has not elapsed, 
                    // then remove a proportion of plant matter
                    // according to the number of time-steps that have elapsed since the spin-up ended

                    int ConstImpactSteps = Convert.ToInt32(humanNPPScenario.Item3 * _Utilities.ConvertTimeUnits("year", globalModelTimestepUnits));

                    if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps)))
                    {
                        //gridCellStocks[actingStock].TotalBiomass -= gridCellStocks[actingStock].TotalBiomass *
                        //    (Math.Min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2)));

                        RemovalRate = (Math.Min(1.0,
                                        Math.Min(((ConstImpactSteps / 12.0) * humanNPPScenario.Item2),
                                        (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2))));
                    }
                }
                else if (humanNPPScenario.Item1 == "temp-escalating-declining")
                {
                    // If the spin-up period has been completed, then apply a level of harvesting
                    // according to the number of time-steps that have elapsed since the spin-up ended
                    if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps)))
                    {
                        //gridCellStocks[actingStock].TotalBiomass -= gridCellStocks[actingStock].TotalBiomass *
                        //    (Math.Min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2)));

                        RemovalRate = Math.Max(0.0, (Math.Min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPScenario.Item2))));
                    }
                    else if ((currentTimestep > (burninSteps + impactSteps)) & (currentTimestep <= (burninSteps + impactSteps + recoverySteps)))
                    {
                        //gridCellStocks[actingStock].TotalBiomass -= gridCellStocks[actingStock].TotalBiomass *
                        //    (Math.Min(1.0, (((burninSteps + impactSteps + recoverySteps - currentTimestep) / 12.0) * humanNPPScenario.Item2)));

                        //RemovalRate = (Math.Min(1.0, (((burninSteps + impactSteps + recoverySteps - currentTimestep) / 12.0) * humanNPPScenario.Item2)));
                        RemovalRate = Math.Max(0.0, Math.Min(1.0, ((int)((impactSteps) - (currentTimestep - (burninSteps + impactSteps))) / 12.0) * humanNPPScenario.Item2));
                    }

                }
                else
                {
                    Debug.Fail("There is no method for the human extraction of NPP scenario specified");
                }

            }

            cellEnvironment["RelativeHANPP"][0] = RemovalRate;
            return(RemovalRate);
        }


        public double FracImpactedHANPPh(double ffor, double fcrop, double furb, double fgra)
        {
            
            double FracImpHANPPh = (fcrop*HANPPh_crop + furb*HANPPh_urban + fgra*HANPPh_grazing)/
                (fcrop * HANPPh_crop + furb * HANPPh_urban + fgra * HANPPh_grazing + ffor * HANPPh_forestry); // weighted averaging to obtain fraction (%)

            //Prevent NaNs arising from division by zero
            if (double.IsNaN(FracImpHANPPh)) FracImpHANPPh = 0.0;

            return FracImpHANPPh;
        }
    }

}
