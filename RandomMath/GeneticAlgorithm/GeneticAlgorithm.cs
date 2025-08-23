using System;
using System.Collections.Generic;
using System.Runtime;
using System.Text;
using System.Threading;

namespace RandomMath.GeneticAlgorithm
{
    /// <summary>
    /// Finds the minimum of some function.  (if you want the max, invert your function)
    /// </summary>
    public class GeneticAlgorithm
    {
        #region Private member variables
        private List<Population> lPopulation;
        private IOptimize function;

        //Setting
        private double MinError;
        private int MaxIt;
        private GASettings settings;

        //Results Tracking
        private Chromosome best;
        private int iterations;
        private double CurrentError;
        private List<double> lAverageFitness;
        private List<double> lBestFitness;
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
        public double BestFitness
        {
            get
            {
                return lPopulation[iterations - 1].BestFitness;
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

        public GeneticAlgorithm(IOptimize inputFunction) : this(inputFunction, new GASettings(100, 0.1, 0.4, 0.8, 1), 0.000001, 100)
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
        public GeneticAlgorithm(IOptimize inputFunction, GASettings settings, double Error, int MaxIterations)
        {
            lBestFitness = new List<double>();
            lAverageFitness = new List<double>();
            lPopulation = new List<Population>();
            iterations = 0;
            this.settings = settings;
            function = inputFunction;
            MinError = Error;
            MaxIt = MaxIterations;

            //Run the first iteration
            Population firstPopulation = new Population(function, settings);

            //Record the 1st iteration, and set Function to the Best
            SaveResults(firstPopulation);
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

            //Pass the previous population to the next generation
            Chromosome[] prevChildren = lPopulation[iterations - 1].CPopulation;
            IOptimize[] functions = lPopulation[iterations - 1].Functions;
            Population newPopulation = new Population(functions, prevChildren, settings);
            
            SaveResults(newPopulation);

            CurrentError = Math.Abs((lPopulation[iterations-1].AverageFitness -
                                     lPopulation[iterations - 2].AverageFitness)
                                    / lPopulation[iterations-1].AverageFitness);

            //Check 
            if ((CurrentError < MinError && CurrentError != 0.0) || iterations > MaxIt)
            {
                PerformNextStep = false;
                iterations--;
            }
            return PerformNextStep;
        }

        private void SaveResults(Population newPopulation)
        {
            //Add the results
            lPopulation.Add(newPopulation);
            lAverageFitness.Add(newPopulation.AverageFitness);
            lBestFitness.Add(newPopulation.BestFitness);
            best = newPopulation.Best;
            iterations++;
        }
        #endregion

        #region private methods


        #endregion
    }

    #region Interface needed for the GeneticAlgorithm

    public interface IOptimize
    {
        int nX { get; }
        
        /// <summary>
        /// Evaluates the function at x, and gives the fitness
        /// </summary>
        /// <param name="X">all inputs will come in at 0-1.  Must scale them out</param>
        /// <returns>the fitness of the function: some numeric value that can be used to compare the function
        /// evaluated at different points.  The fitness will be minimized</returns>
        double Eval(double[] X);

        IOptimize DeepCopy();
    }

    #endregion

    public class GASettings
    {
        public int NumSamples { get; set; }
        public double MutationProbability { get; set; }
        public double PurturbationProbability { get; set; }
        public double CrossoverProbability { get; set; }
        public int NumCores { get; set; }

        public GASettings(int numSamples, double mutationProbability, double crossoverProbability, double purturbationProbability, int numCores)
        {
            //Make sure nsamples is even
            NumSamples = ((int)(numSamples / 2.0d)) * 2;
            MutationProbability = mutationProbability;
            CrossoverProbability = crossoverProbability;
            NumCores = numCores;
            PurturbationProbability = purturbationProbability;
        }
    }
    public static class ThreadRandom
    {
        private static readonly ThreadLocal<Random> threadRandom =
            new ThreadLocal<Random>(() => new Random(Guid.NewGuid().GetHashCode()));

        public static Random GetThreadRandom()
        {
            return threadRandom.Value;
        }
        public static double NextDouble()
        {
            return GetThreadRandom().NextDouble();
        }

        public static int Next(int maxValue)
        {
            return GetThreadRandom().Next(maxValue);
        }
    }
}