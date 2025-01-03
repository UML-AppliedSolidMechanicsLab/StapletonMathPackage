using System;
using System.Collections.Generic;
using System.Text;

namespace RandomMath.GeneticAlgorithm
{
    public class Population
    {
        #region private members
        private Chromosome[] cPopulation;
        private Chromosome[] cChildren;
        private IOptimize function;
        private int nSamples;
        private double probMutate;
        private double probXOver;
        private Random random = new Random();
        private double AveFitness;
        #endregion

        #region Public Members
        public Chromosome Best
        {
            get
            {
                return cPopulation[0];
            }
        }
        public Chromosome[] Children
        {
            get
            {
                return cChildren;
            }
        }
        public double AverageFitness
        {
            get
            {
                return AveFitness;
            }
        }
        #endregion

        #region Constructors
        public Population(IOptimize Function, int NumberOfSamples, double ProbabilityOfMutation, double ProbabilityOfXOver)
        {

            //Make sure nSamples is even
            nSamples = ((int)(NumberOfSamples / 2d)) * 2;
            Chromosome[] firstPopulation = FirstPopulation(Function, nSamples, ProbabilityOfMutation);
            Initiate(Function, firstPopulation, NumberOfSamples, ProbabilityOfMutation, ProbabilityOfXOver);
        }

        public Population(IOptimize Function, Chromosome[] PrevChildren, int NumberOfSamples, double ProbabilityOfMutation, double ProbabilityOfXOver)
        {
            Initiate(Function, PrevChildren, NumberOfSamples, ProbabilityOfMutation, ProbabilityOfXOver);
        }

        private void Initiate(IOptimize Function, Chromosome[] PrevChildren, int NumberOfSamples, double ProbabilityOfMutation, double ProbabilityOfXOver)
        {
            function = Function;
            cPopulation = PrevChildren;
            probMutate = ProbabilityOfMutation;
            probXOver = ProbabilityOfXOver;
            //Make sure nSamples is even
            nSamples = ((int)(NumberOfSamples / 2.0d)) * 2;

            cChildren = new Chromosome[nSamples];

            FindAverageFitness();

            MatingTournament();

        }

        #endregion

        #region private methods
        private Chromosome[] FirstPopulation(IOptimize Function, int NumberOfSamples, double propMutate)
        {
            Chromosome[] firstPop = new Chromosome[NumberOfSamples];

            for (int i = 0; i < NumberOfSamples; i++)
            {

                double[] inputString = new double[Function.nX];

                for (int j = 0; j < Function.nX; j++)
                {

                    inputString[j] = random.NextDouble();

                }
                firstPop[i] = new Chromosome(Function, inputString, propMutate);
            }
            return firstPop;
        }

        private void MatingTournament()
        {

            for (int i = 0; i < nSamples / 2; i++)
            {

                //Chose 4 samples for tournament
                int[] randomIndices = new int[4];
                for (int j = 0; j < 4; j++)
                {
                    randomIndices[j] = RandomIndex(nSamples);
                }
                //Since they are all sorted, the lowest 2 #'s are the parents
                Array.Sort(randomIndices);

                cChildren[2 * i] = new Chromosome(cPopulation[randomIndices[0]]);
                cChildren[2 * i + 1] = new Chromosome(cPopulation[randomIndices[1]]);

                //XOver and Mutate
                BabyMaker(cChildren[2 * i], cChildren[2 * i + 1]);
            }
        }

        private void BabyMaker(Chromosome child1, Chromosome child2)
        {

            bool firstHasBeenChanged = false;
            bool secondHasBeenChanged = false;

            //perform xover twice
            if (random.NextDouble() < probXOver)
            {

                firstHasBeenChanged = true;
                secondHasBeenChanged = true;

                XOver(child1.String, child2.String);
                XOver(child1.String, child2.String);

            }
            //Mutate
            child1.Mutate(firstHasBeenChanged);
            child2.Mutate(secondHasBeenChanged);

        }

        private int RandomIndex(int nPoints)
        {

            double rN = random.NextDouble();
            double rIndex = (nPoints * rN);
            int rInt = (int)(rIndex);
            return rInt;

        }

        private void XOver(double[] c1, double[] c2)
        {

            double[] tempC1 = new double[c1.Length];
            int ranIndex = RandomIndex(c1.Length);

            for (int i = 0; i < c1.Length; i++)
            {
                tempC1[i] = c1[i];
            }
            for (int i = ranIndex; i < c1.Length; i++)
            {
                c1[i] = c2[i];// + (c1[i] - c2[i])/2.0*random.NextDouble();
                c2[i] = tempC1[i];// + (tempC1[i] - c2[i])/2.0*random.NextDouble();;
            }

        }

        private void FindAverageFitness()
        {

            Array.Sort(cPopulation);

            AveFitness = cPopulation[0].Fitness;

            for (int i = 1; i < nSamples; i++)
            {

                AveFitness += cPopulation[i].Fitness;
            }

            AveFitness = AveFitness / nSamples;
        }


        #endregion

        #region Public Methods

        public Chromosome[] NewPopulation()
        {

            //Combine children and population, and keep nsamples of the best ones
            Chromosome[] BigHappyFamily = new Chromosome[2 * nSamples];

            //copy manually to ensure deep copy!!
            cPopulation.CopyTo(BigHappyFamily, 0);
            cChildren.CopyTo(BigHappyFamily, nSamples);

            //now, sort the big happy family
            Array.Sort(BigHappyFamily);
            //now take the better half
            Chromosome[] newPop = new Chromosome[nSamples];
            for (int i = 0; i < nSamples; i++)
            {

                newPop[i] = new Chromosome(BigHappyFamily[i]);
            }

            return newPop;
        }

        #endregion
    }



    /// <summary>
    /// Description of Chromosome.
    /// </summary>
    public class Chromosome : IComparable
    {
        #region private members
        private IOptimize function;
        private double[] dString;
        private double fitness;
        private double pMutate;
        private Random random = new Random();
        #endregion

        #region public members
        public double[] String
        {
            get
            {
                return dString;
            }
            set
            {
                dString = value;
            }
        }
        public double Fitness
        {
            get
            {
                return fitness;
            }
        }
        #endregion

        #region constructor
        public Chromosome(IOptimize Function, double[] inputString, double ProbabilityOfMutation)
        {
            function = Function;
            fitness = function.Eval(inputString);
            pMutate = ProbabilityOfMutation;
            dString = inputString;

        }
        //Put this constructor to copy a Chromosome
        public Chromosome(Chromosome ChromosomeToCopy)
        {
            function = ChromosomeToCopy.function;
            fitness = ChromosomeToCopy.fitness;
            pMutate = ChromosomeToCopy.pMutate;
            dString = new double[ChromosomeToCopy.String.Length];
            for (int i = 0; i < ChromosomeToCopy.String.Length; i++)
            {
                dString[i] = ChromosomeToCopy.dString[i];
            }

        }
        #endregion

        #region public methods
        public int CompareTo(object o)
        {
            Chromosome c = (Chromosome)o;
            //This will spit out a negative interger if this is "smaller".  Since smaller is good, smaller numerically is "smaller" in the icomparable sense
            //return (int)(c.fitness - fitness);
            double dReturn = fitness - c.fitness;
            int iReturn = 0;
            if (dReturn > 0)
            {
                iReturn = 1;
            }
            else if (dReturn < 0)
            {
                iReturn = -1;
            }

            return iReturn;

        }

        public void Mutate(bool hasBeenChanged)
        {

            //Fine tuning Mutation		
            for (int i = 0; i < function.nX; i++)
            {
                double r = random.NextDouble();
                if (r < pMutate)
                {

                    hasBeenChanged = true;

                    if ((1.0 - dString[i]) < .005)
                    {

                        dString[i] -= 0.005;

                    }
                    else if ((dString[i]) < .005)
                    {

                        dString[i] += 0.005;

                    }
                    else
                    {
                        r = random.NextDouble();
                        if (r < 0.5)
                        {
                            dString[i] += 0.005;
                        }
                        else dString[i] -= 0.005;

                    }
                }
            }
            //A second type of mutation
            for (int i = 0; i < function.nX; i++)
            {

                double r = random.NextDouble();

                if (r < pMutate)
                {

                    hasBeenChanged = true;

                    dString[i] = random.NextDouble();

                }
            }

            //now evaluate it again if it has been changed
            if (hasBeenChanged)
            {

                fitness = function.Eval(dString);
            }
        }
        #endregion
    }
}