using System;
using System.Collections.Generic;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using static RandomMath.RootFinding;

namespace RandomMath.GeneticAlgorithm
{
    public class Population
    {
        #region private members
        private Chromosome[] cPopulation;
        private double aveFitness;
        private double bestFitness;
        private IOptimize[] functions;

        #endregion

        #region Public Members
        public Chromosome Best
        {
            get
            {
                return cPopulation[0];
            }
        }

        public double AverageFitness
        {
            get
            {
                return aveFitness;
            }
        }

        public double BestFitness
        {
            get
            {
                return bestFitness;
            }
        }

        public Chromosome[] CPopulation
        {
            get
            {
                return cPopulation;
            }
        }

        public IOptimize[] Functions
        {
            get
            {
                return functions;
            }
        }
        #endregion

        #region Constructors
        /// <summary>
        /// Constructor for the first population only
        /// </summary>
        /// <param name="function">function to be evaluated</param>
        public Population(IOptimize function, GASettings settings)
        {
            Chromosome[] firstPopulation = MakeFirstPopulation(function, settings);

            //Make a function for each core so that there is no conflict if running on multiple cores
            functions = new IOptimize[settings.NumCores];
            for (int i = 0; i < settings.NumCores; i++)
            {
                functions[i] = function.DeepCopy();
            }

            //Evaluate all of the functions
            //Now do the evaluation part in paralell for longer evaluations
            Parallel.For(0, firstPopulation.Length, new ParallelOptions { MaxDegreeOfParallelism = settings.NumCores }, i =>
            {
                int coreIndex = i % settings.NumCores; // Assign optimizer per thread (round-robin)
                firstPopulation[i].MutateAndEvaluateFitness(true, 0.0, 0.0, functions[coreIndex]);
            });

            //Save the results
            cPopulation = firstPopulation;
            Array.Sort(cPopulation); // Sort the population based on fitness
            FindAverageFitness();
            bestFitness = cPopulation[0].Fitness;
        }

        public Population(IOptimize[] functions, Chromosome[] PrevChildren, GASettings settings)
        {
            //Save previous children as parents

            this.functions = functions;

            Chromosome[] cChildren = MatingTournaments(settings, PrevChildren);

            cPopulation = CombineBestOfChildrenAndParentsAndSort(cChildren, PrevChildren);

            FindAverageFitness();

            bestFitness = cPopulation[0].Fitness;
        }

        #endregion

        #region Public Methods

        #endregion

        #region private methods
        private Chromosome[] MakeFirstPopulation(IOptimize function, GASettings settings)
        {
            Chromosome[] firstPop = new Chromosome[settings.NumSamples];

            for (int i = 0; i < settings.NumSamples; i++)
            {

                double[] inputString = new double[function.nX];

                for (int j = 0; j < function.nX; j++)
                {

                    inputString[j] = ThreadRandom.NextDouble();

                }
                firstPop[i] = new Chromosome(inputString);
            }
            return firstPop;
        }

        private Chromosome[] MatingTournaments(GASettings settings, Chromosome[] cParents)
        {
            int nPairs = cParents.Length / 2;
            int nCores = settings.NumCores;

            Chromosome[] cChildren = new Chromosome[settings.NumSamples];

            //First, make the children population so that my parallel loop is threadsafe:
            for (int i = 0; i < nPairs; i++)
            {
                //Chose 4 samples for tournament
                int[] randomIndices = new int[4];
                for (int j = 0; j < 4; j++)
                {
                    randomIndices[j] = ThreadRandom.Next(cParents.Length);
                }
                //Since they are all sorted, the lowest 2 #'s are the parents
                Array.Sort(randomIndices);

                cChildren[2 * i] = cParents[randomIndices[0]].DeepCopy();
                cChildren[2 * i + 1] = cParents[randomIndices[1]].DeepCopy();
            }

            //Now do the evaluation part in paralell for longer evaluations
            Parallel.For(0, nPairs, new ParallelOptions { MaxDegreeOfParallelism = nCores }, i =>
            {
                int coreIndex = i % nCores; // Assign optimizer per thread (round-robin)

                //XOver and MutateAndEvaluateFitness
                BabyMaker(functions[coreIndex], settings.CrossoverProbability, settings.MutationProbability,
                    settings.PurturbationProbability, cChildren[2 * i], cChildren[2 * i + 1]);

            });

            return cChildren;
        }

        private void BabyMaker(IOptimize function, double probXOver, double probMutation,
            double probPurturbation, Chromosome child1, Chromosome child2)
        {

            bool haveBeenChanged = false;

            //perform xover twice
            if (ThreadRandom.NextDouble() < probXOver)
            {
                haveBeenChanged = true;

                XOver(child1.String, child2.String);
                //XOver(child1.String, child2.String);
            }
            //MutateAndEvaluateFitness
            child1.MutateAndEvaluateFitness(haveBeenChanged, probMutation, probMutation, function);
            child2.MutateAndEvaluateFitness(haveBeenChanged, probMutation, probMutation, function);

        }

        private void XOver(double[] c1, double[] c2)
        {

            double[] tempC1 = new double[c1.Length];
            Array.Copy(c1, tempC1, c1.Length);

            int ranIndex = ThreadRandom.Next(c1.Length);

            for (int i = ranIndex; i < c1.Length; i++)
            {
                c1[i] = c2[i];// + (c1[i] - c2[i])/2.0*random.NextDouble();
                c2[i] = tempC1[i];// + (tempC1[i] - c2[i])/2.0*random.NextDouble();;
            }

        }

        /// <summary>
        /// Creates a new population from the current population and the children, 
        /// taking the best half of both
        /// </summary>
        /// <returns>population of the best half of parents and children, sorted</returns>
        private Chromosome[] CombineBestOfChildrenAndParentsAndSort(Chromosome[] cChildren, Chromosome[] cParents)
        {

            //Combine children and population, and keep nsamples of the best ones
            Chromosome[] parentsAndChildren = new Chromosome[2 * cChildren.Length];

            //copy manually to ensure deep copy!!
            cParents.CopyTo(parentsAndChildren, 0);
            cChildren.CopyTo(parentsAndChildren, cParents.Length);

            //now, sort the big happy family
            Array.Sort(parentsAndChildren);
            //now take the "fitter" half
            Chromosome[] newPop = new Chromosome[cChildren.Length];

            for (int i = 0; i < cChildren.Length; i++)
            {
                newPop[i] = parentsAndChildren[i].DeepCopy();
            }

            return newPop;
        }

        private void FindAverageFitness()
        {
            double sum = 0;

            for (int i = 0; i < cPopulation.Length; i++)
            {
                sum += cPopulation[i].Fitness; ;
            }
            aveFitness = sum / cPopulation.Length;
        }

        #endregion

    }
}