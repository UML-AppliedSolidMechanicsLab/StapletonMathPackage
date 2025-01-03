using System;
using System.Collections.Generic;
using System.Text;

namespace RandomMath.GeneticAlgorithm
{
    /// <summary>
    /// Description of GeneticAlgorithm.  Finds the minimum of some function.  (if you want the max, invert your function)
    /// </summary>
    public class GeneticAlgorithm
    {
        #region Private member variables
        private List<Population> lPopulation;
        private Chromosome best;
        private int iterations;
        private double CurrentError;
        private List<double> lAverageFitness;
        private List<double> lBestFitness;
        private IOptimize function;
        private int nSamples;
        private double probMutate;
        private double probXOver;
        private double MinError;
        private int MaxIt;
        #endregion

        #region Public Properties
        public int Iterations
        {
            get
            {
                return iterations;
            }
        }
        public Chromosome Best
        {
            get
            {
                return best;
            }
        }
        public double RelError
        {
            get
            {
                return CurrentError;
            }
        }
        public List<double> AverageFitnesses
        {
            get
            {
                return lAverageFitness;
            }
        }
        public List<double> BestFitnesses
        {
            get
            {
                return lBestFitness;
            }
        }
        public IOptimize Function
        {
            get
            {
                return function;
            }
        }
        public double MaximumIterations
        {
            get
            {
                return MaxIt;
            }
        }
        #endregion

        #region Constructors

        public GeneticAlgorithm(IOptimize inputFunction) : this(inputFunction, 100, 0.1, 0.8, 0.000001, 100)
        {

        }

        #region
        /// <summary>
        /// GeneticAlgorithm uses a genetic algorithm (duh) to optimize some function of interface IOptimize.
        /// </summary>
        /// <param name="Function">function to be optimized, along with method to compare 2 function evaluations</param>
        /// <param name="nSamples">number of chromosome samples per population</param>
        /// <param name="probMutate">the probability that a string in a chromosome will mutate (0 - 1) </param>
        /// <param name="probXOver">the probability that a crossover between children will happen (0 - 1) </param>
        /// <param name="Error">The relative error between best chromosome, which dictates whether optimal is found</param>
        #endregion
        public GeneticAlgorithm(IOptimize inputFunction, int NumberOfSamples, double ProbMutate, double ProbXOver, double Error, int MaxIterations)
        {
            lBestFitness = new List<double>();
            lAverageFitness = new List<double>();
            lPopulation = new List<Population>();
            iterations = 0;
            nSamples = NumberOfSamples;
            function = inputFunction;
            probMutate = ProbMutate;
            probXOver = ProbXOver;
            MinError = Error;
            MaxIt = MaxIterations;

            //Run the first iteration
            lPopulation.Add(new Population(function, nSamples, probMutate, probXOver));

            //Record the 1st iteration, and set Function to the Best
            lAverageFitness.Add(lPopulation[iterations].AverageFitness);
            lBestFitness.Add(function.Eval(lPopulation[iterations].Best.String));
            best = lPopulation[iterations].Best;

            iterations++;
        }

        #endregion

        #region Private Methods


        #endregion

        #region Public Methods

        public void Run()
        {
            bool PerformNextStep = true;


            while (PerformNextStep)
            {

                PerformNextStep = RunSingleStep();

            }
        }

        //This version is for use with a bagroundworker thread
        public bool RunSingleStep()
        {
            bool PerformNextStep = true;

            Chromosome[] newPopulation = lPopulation[iterations - 1].NewPopulation();
            lPopulation.Add(new Population(function, newPopulation, nSamples, probMutate, probXOver));

            CurrentError = Math.Abs((lPopulation[iterations].AverageFitness -
                                     lPopulation[iterations - 1].AverageFitness)
                                    / lPopulation[iterations].AverageFitness);
            //Add the results
            lAverageFitness.Add(lPopulation[iterations].AverageFitness);
            lBestFitness.Add(function.Eval(lPopulation[iterations].Best.String));
            best = lPopulation[iterations].Best;
            iterations++;

            //Check 
            if ((CurrentError < MinError && CurrentError != 0.0) || iterations > MaxIt)
            {

                PerformNextStep = false;
                iterations--;
            }
            return PerformNextStep;
        }
        #endregion

    }

    #region Interface needed for the GeneticAlgorithm

    public interface IOptimize
    {
        int nX { get; }
        double[] currentC { get; }
        /// <summary>
        /// Evaluates the function at x, and gives the fitness
        /// </summary>
        /// <param name="X">all inputs will come in at 0-1.  Must scale them out</param>
        /// <returns>the fitness of the function: some numeric value that can be used to compare the function
        /// evaluated at different points</returns>
        double Eval(double[] X);

    }

    #endregion
}