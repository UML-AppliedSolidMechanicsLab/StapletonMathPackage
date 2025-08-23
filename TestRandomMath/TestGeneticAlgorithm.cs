using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using RandomMath.GeneticAlgorithm;
using NUnit.Framework.Internal;


namespace TestRandomMath
{
    [TestFixture]
    public class TestGeneticAlgorithm
    {
        [Test]
        public void GA_Minimize_1Variables()
        {
            // Define the function to be optimized
            IOptimize f = new FunToMin_1Var();
            double minError = 0.0000001;
            int maxIterations = 1000;
            GASettings settings = new GASettings(500, 0.01, 0.0, 0.5, 1);
            GeneticAlgorithm ga = new GeneticAlgorithm(f, settings, minError, maxIterations);
            ga.Run();
            double y = ga.Best.Fitness;
            double x = Math.Abs(ga.Best.String[0] - 0.7);

            Assert.IsNotNull(x, "Best chromosome should not be null.");
            Assert.Less(y, 0.0005, "Fitness should be close to 0 for a simple quadratic.");
            Assert.Less(x, 0.0005, "X should be close to 0.7 for this one.");
            Assert.That(x, Is.EqualTo(y).Within(1e-6), "Fitness and |x[0] - 0.7| should match exactly.");
        }

        [Test]
        public void GA_Minimizes_5Variables()
        {
            // Arrange
            var function = new FunToMin_5Var();
            var settings = new GASettings(numSamples: 500, mutationProbability: 0.4, crossoverProbability: 0.8, 
                purturbationProbability: 0.8, numCores: 4);
            var ga = new GeneticAlgorithm(function, settings, Error: 1e-8, MaxIterations: 200);

            // Act
            ga.Run();

            // Assert
            var best = ga.Best;
            Assert.IsNotNull(best, "Best chromosome should not be null.");
            Assert.Less(best.Fitness, 0.0001, "Fitness should be close to 0 for a simple quadratic.");

            for (int i = 1; i < ga.BestFitnesses.Count; i++)
            {
                Assert.LessOrEqual(ga.BestFitnesses[i], ga.BestFitnesses[i - 1], $"Fitness increased at generation {i}");
            }
        }
    }
    public class FunToMin_1Var : IOptimize
    {
        public int nX => 1; // 1-dimensional input
        public double Eval(double[] x)
        {
            // Example function: f(x) = (x-0.7): minimum is 0.7
            return Math.Abs(x[0] - 0.7);
        }
        public IOptimize DeepCopy()
        {
            return new FunToMin_1Var();
        }

    }

    public class FunToMin_5Var : IOptimize
    {
        public int nX => 5; // 5-dimensional input

        public double Eval(double[] X)
        {
            double sum = 0;
            for (int i = 0; i < X.Length; i++)
            {
                sum += X[i] * X[i];
            }
            return sum; // lower is better
        }

        public IOptimize DeepCopy()
        {
            return new FunToMin_5Var();
        }
    }
}
