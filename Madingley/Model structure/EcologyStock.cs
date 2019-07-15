using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

using System.IO;

namespace Madingley
{
    /// <summary>
    /// A class to specify, initialise and run ecological processes pertaining to stocks
    /// </summary>
    class EcologyStock
    {
        /// <summary>
        /// An instance of the Autotroph Processor for this model
        /// </summary>
        AutotrophProcessor MarineNPPtoAutotrophStock;

        /// <summary>
        /// An instance of the plant model class
        /// </summary>
        RevisedTerrestrialPlantModel DynamicPlantModel;

        /// <summary>
        /// An instance of the class for human appropriation of NPP
        /// </summary>
        HumanAutotrophMatterAppropriation HANPP;



        public void InitializeEcology()
        {
            //Initialize the autotrophprocessor
            MarineNPPtoAutotrophStock = new AutotrophProcessor();

            // Initialise the plant model
            DynamicPlantModel = new RevisedTerrestrialPlantModel();

            // Initialise the human NPP appropriation class
            HANPP = new HumanAutotrophMatterAppropriation();

        }


        /// <summary>
        /// Run ecological processes that operate on stocks within a single grid cell
        /// </summary>
        ///<param name="gridCellStocks">The stocks in the current grid cell</param>
        ///<param name="actingStock">The acting stock</param>
        ///<param name="cellEnvironment">The stocks in the current grid cell</param>
        ///<param name="environmentalDataUnits">List of units associated with the environmental variables</param>
        ///<param name="humanNPPScenario">The human appropriation of NPP scenario to apply</param>
        ///<param name="madingleyStockDefinitions">The functional group definitions for stocks in the model</param>
        ///<param name="currentTimeStep">The current model time step</param>
        ///<param name="burninSteps">The number of time steps to spin the model up for before applying human impacts</param>
        ///<param name="impactSteps">The number of time steps to apply human impacts for</param>
        ///<param name="globalModelTimeStepUnit">The time step unit used in the model</param>
        ///<param name="trackProcesses">Whether to track properties of ecological processes</param>
        ///<param name="tracker">An instance of the ecological process tracker</param>
        ///<param name="globalTracker">An instance of the global process tracker</param>
        ///<param name="currentMonth">The current model month</param>
        ///<param name="outputDetail">The level of detail to use in outputs</param>
        ///<param name="specificLocations">Whether to run the model for specific locations</param>
        ///<param name="impactCell">Whether this cell should have human impacts applied</param>
        public void RunWithinCellEcology(GridCellStockHandler gridCellStocks, int[] actingStock, SortedList<string, double[]> cellEnvironment,
            SortedList<string, string> environmentalDataUnits, Tuple<string, double, double> humanNPPScenario, 
            FunctionalGroupDefinitions madingleyStockDefinitions, 
            uint currentTimeStep, uint burninSteps, uint impactSteps,uint recoverySteps, uint instantStep, uint numInstantSteps, string globalModelTimeStepUnit, Boolean trackProcesses, 
            ProcessTracker tracker, 
            GlobalProcessTracker globalTracker, uint currentMonth, 
            string outputDetail, bool specificLocations, Boolean impactCell)
        {
            double loss = 0.0;
            double gain = 0.0;

            int ScenarioYear;
            if(currentTimeStep < burninSteps)
            {
                ScenarioYear = 0;
                if (madingleyStockDefinitions.GetTraitNames("impact state", actingStock[0]) == "primary")
                {
                    gridCellStocks[actingStock].FractionalArea = cellEnvironment["Fprimary"][ScenarioYear];
                }
                else if (madingleyStockDefinitions.GetTraitNames("impact state", actingStock[0]) == "secondary")
                {
                    gridCellStocks[actingStock].FractionalArea = cellEnvironment["Fsecondary"][ScenarioYear];
                }
                else
                {
                    //All HANPPlc comes from impacted lands
                    gridCellStocks[actingStock].FractionalArea = 1 - (cellEnvironment["Fprimary"][ScenarioYear] + cellEnvironment["Fsecondary"][ScenarioYear]);
                }
            }
            else
            {
                ScenarioYear = (int)Math.Floor((currentTimeStep-burninSteps) / 12.0);



				if (humanNPPScenario.Item1 == "ssp") // annual hanpp (gC/m2/yr) input
				{

					//NEEDS amending as is only set up for terrestrial cells at present.
					//if (madingleyStockDefinitions.GetTraitNames("Realm", actingStock[0]) == "marine")
					//{
					//    // Run the autotroph processor
					//    MarineNPPtoAutotrophStock.ConvertNPPToAutotroph(cellEnvironment, gridCellStocks, actingStock, environmentalDataUnits["LandNPP"], 
					//        environmentalDataUnits["OceanNPP"], currentTimeStep,globalModelTimeStepUnit,tracker,globalTracker ,outputDetail,specificLocations,currentMonth);
					//}
					//else if (madingleyStockDefinitions.GetTraitNames("Realm", actingStock[0]) == "terrestrial")


					if (madingleyStockDefinitions.GetTraitNames("impact state", actingStock[0]) == "primary")
					{
						loss = cellEnvironment["Primary loss"][ScenarioYear] / 12.0;
						//calculate the change in total biomass as a result of coverage changes
						// assumes that the biomass density stays the same as coverage goes down (ie if chopping down some forest - the density of the remaining stays the same)
						// However if coverage increases, the density declines, as biomass is only added by NPP..
						//if (cellEnvironment["Primary loss"][ScenarioYear] < gridCellStocks[actingStock].FractionalArea)

					}
					else if (madingleyStockDefinitions.GetTraitNames("impact state", actingStock[0]) == "secondary")
					{
						loss = cellEnvironment["Secondary loss"][ScenarioYear] / 12.0;
						gain = cellEnvironment["Secondary gain"][ScenarioYear] / 12.0;

					}
					else if (madingleyStockDefinitions.GetTraitNames("impact state", actingStock[0]) == "impacted")
					{

						loss = cellEnvironment["Secondary gain"][ScenarioYear] / 12.0;
						gain = (cellEnvironment["Secondary loss"][ScenarioYear] + cellEnvironment["Primary loss"][ScenarioYear]) / 12.0;

					}

					if (gridCellStocks[actingStock].FractionalArea.CompareTo(0.0) > 0)
						gridCellStocks[actingStock].TotalBiomass *= 1.0 - Math.Min(1.0, loss / gridCellStocks[actingStock].FractionalArea);
					gridCellStocks[actingStock].FractionalArea = Math.Max(0.0, gridCellStocks[actingStock].FractionalArea - loss + gain);

				}


				// CRAFTY HANPP 
				else if (humanNPPScenario.Item1 == "hanpp_crafty") // annual hanpp (gC/m2/yr) input

				{
					// @todo consider impact state.. 

					if (madingleyStockDefinitions.GetTraitNames("Realm", actingStock[0]) == "marine")
					{
						// Run the autotroph processor
						MarineNPPtoAutotrophStock.ConvertNPPToAutotroph(cellEnvironment, gridCellStocks, actingStock, environmentalDataUnits["LandNPP"],
							environmentalDataUnits["OceanNPP"], currentTimeStep, globalModelTimeStepUnit, tracker, globalTracker, outputDetail, specificLocations, currentMonth);
					}
					else if (madingleyStockDefinitions.GetTraitNames("Realm", actingStock[0]) == "terrestrial")
					{

		 				// have to have another NC file capsuling aft_fractions (17 types) at 0.25 deg grid. 
						double standingstock_remained_perc = cellEnvironment["remained_standingstock_perc"][ScenarioYear]; // remaind standing stock (%/yr)
						double monthly_loss = (1 - standingstock_remained_perc) / 12.0; 
			
                        // Apply standing stock changes
                        if (monthly_loss > 0)
							gridCellStocks[actingStock].TotalBiomass *= 1.0 - Math.Min(1.0, monthly_loss);

                        // @todo fractional area change

					}
				}
            }

            // Run the dynamic plant model to update the leaf stock for this time step
            double WetMatterNPP = DynamicPlantModel.UpdateLeafStock(cellEnvironment, gridCellStocks, actingStock, currentTimeStep, madingleyStockDefinitions.
                GetTraitNames("leaf strategy", actingStock[0]).Equals("deciduous"), globalModelTimeStepUnit, tracker, globalTracker, currentMonth,
                outputDetail, specificLocations);
                        
            double fhanpp = HANPP.RemoveHumanAppropriatedMatter(WetMatterNPP, cellEnvironment, humanNPPScenario, gridCellStocks, actingStock,
                currentTimeStep, ScenarioYear, burninSteps, impactSteps, recoverySteps, instantStep, numInstantSteps, impactCell, globalModelTimeStepUnit, madingleyStockDefinitions,
                DynamicPlantModel.CalculateFracEvergreen(cellEnvironment["Fraction Year Frost"][currentTimeStep]));



            // Apply human appropriation of NPP
            gridCellStocks[actingStock].TotalBiomass += WetMatterNPP * (1.0 - fhanpp);
            if (globalTracker.TrackProcesses)
            {
                globalTracker.RecordHANPP((uint)cellEnvironment["LatIndex"][0], (uint)cellEnvironment["LonIndex"][0], (uint)actingStock[0],
                    fhanpp);
            }

            if (gridCellStocks[actingStock].TotalBiomass < 0.0) gridCellStocks[actingStock].TotalBiomass = 0.0;

            //}
            //else
            //{
            //    Debug.Fail("Stock must be classified as belonging to either the marine or terrestrial realm");
            //}
        }
    }
}
