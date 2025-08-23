using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RandomMath.GeneticAlgorithm
{
    /// <summary>
    /// Basic building block of a design.  Carries the dString, which is the 
    /// input variables (x) and the fitness, which is how that design fared
    /// in the evaluation.
    /// </summary>
    public class Chromosome : IComparable
    {
        #region private members
        private double[] dString;
        private double fitness;
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
        public Chromosome(double[] inputString)
        {
            dString = inputString;
        }

        #endregion

        #region public methods

        public Chromosome DeepCopy()
        {
            double[] copiedString = new double[dString.Length];
            Array.Copy(dString, copiedString, dString.Length);

            Chromosome copiedChromosome = new Chromosome(copiedString);
            copiedChromosome.fitness = fitness;
            return copiedChromosome;
        }
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
        public void MutateAndEvaluateFitness(bool hasBeenChanged, double probabilityOfMutation, double probOfPurturbation, IOptimize function)
        {	
            for (int i = 0; i < dString.Length; i++)
            {
                double r = ThreadRandom.GetThreadRandom().NextDouble();
                if (r < probabilityOfMutation)
                {
                    hasBeenChanged = true;
                    // Large mutation: completely new value between 0 and 1
                    dString[i] = ThreadRandom.NextDouble();
                }

                r = ThreadRandom.GetThreadRandom().NextDouble();
                if (r < probOfPurturbation)
                {
                    hasBeenChanged = true;
                    // Large mutation: completely new value between 0 and 1
                    // Small mutation: add or subtract up to 0.001
                    double delta = (ThreadRandom.NextDouble() * 0.001);


                    if (ThreadRandom.NextDouble() < 0.5)
                    { delta *= -1; }

                    double value = dString[i] + delta;
                    // Ensure value is between 0 and 1
                    dString[i] = Math.Max(0.0, Math.Min(1.0, value));
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
