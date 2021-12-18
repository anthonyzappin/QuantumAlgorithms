
namespace Microsoft.Quantum.Samples.IntegerFactorization {
    open Microsoft.Quantum.Random;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Diagnostics;

    operation FactorSemiprimeInteger(number : Int, useRobustPhaseEstimation : Bool) 
    : (Int, Int) {
        // First check the most trivial case, if the provided number is even
        if (number % 2 == 0) {
            Message("An even number has been given; 2 is a factor.");
            return (number / 2, 2);
        }

        mutable foundFactors = false;
        mutable factors = (1, 1);

        repeat {

            let generator = DrawRandomInt(1, number - 1);

            if (IsCoprimeI(generator, number)) {

                Message($"Estimating period of {generator}");

                let period = EstimatePeriod(generator, number, useRobustPhaseEstimation);
                set (foundFactors, factors) = MaybeFactorsFromPeriod(number, generator, period);
            }
            else {
                let gcd = GreatestCommonDivisorI(number, generator);
                
                Message($"We have guessed a divisor of {number} to be {gcd} by accident.");

                set foundFactors = true;
                set factors = (gcd, number / gcd);
            }
        } 
        until foundFactors
        fixup {
            Message("The estimated period did not yield a valid factor, trying again.");
        }
        
        return factors;
    }

    operation ApplyOrderFindingOracle(
        generator : Int, modulus : Int, power : Int, target : Qubit[]
    )
    : Unit
    is Adj + Ctl {
        Fact(IsCoprimeI(generator, modulus), "`generator` and `modulus` must be co-prime");
        MultiplyByModularInteger(ExpModI(generator, power, modulus), modulus, LittleEndian(target));
    }

    operation EstimatePeriod(
        generator : Int, modulus : Int, useRobustPhaseEstimation : Bool
    ) 
    : Int {
        Fact(IsCoprimeI(generator, modulus), "`generator` and `modulus` must be co-prime");

        mutable result = 1;
        let bitsize = BitSizeI(modulus);
        let bitsPrecision = 2 * bitsize + 1;

        mutable frequencyEstimate = 0;

        repeat {

            set frequencyEstimate = EstimateFrequency(
                generator, modulus, useRobustPhaseEstimation, bitsize 
            );

            if (frequencyEstimate != 0) {
                set result = PeriodFromFrequency(modulus,frequencyEstimate, bitsPrecision, result);
            }
            else {
                Message("The estimated frequency was 0, trying again.");
            }
        }
        until(ExpModI(generator, result, modulus) == 1)
        fixup {
            Message("The estimated period from continued fractions failed, trying again.");
        }

        return result;
    }

    operation EstimateFrequency(
        generator : Int, 
        modulus : Int,
        useRobustPhaseEstimation : Bool, 
        bitsize : Int
    )
    : Int {
        mutable frequencyEstimate = 0;
        let bitsPrecision =  2 * bitsize + 1;
        
        use eigenstateRegister = Qubit[bitsize];

        let eigenstateRegisterLE = LittleEndian(eigenstateRegister);
        ApplyXorInPlace(1, eigenstateRegisterLE);

        let oracle = DiscreteOracle(ApplyOrderFindingOracle(generator, modulus, _, _));

        if (useRobustPhaseEstimation) {

            let phase = RobustPhaseEstimation(bitsPrecision, oracle, eigenstateRegisterLE!);
            
            set frequencyEstimate = Round(((phase * IntAsDouble(2 ^ bitsPrecision)) / 2.0) / PI());
        }
        else {
            use register = Qubit[bitsPrecision];
            let frequencyEstimateNumerator = LittleEndian(register);  

            QuantumPhaseEstimation(
                oracle, eigenstateRegisterLE!, LittleEndianAsBigEndian(frequencyEstimateNumerator)
            ); 
            
            set frequencyEstimate = MeasureInteger(frequencyEstimateNumerator);
        }
        
        ResetAll(eigenstateRegister);

        return frequencyEstimate;
    }

    function PeriodFromFrequency(
        modulus : Int, 
        frequencyEstimate : Int, 
        bitsPrecision : Int, 
        currentDivisor : Int
    )
    : Int {
        
        let (numerator, period) = (ContinuedFractionConvergentI(Fraction(frequencyEstimate, 2 ^ bitsPrecision), modulus))!;
        
        let (numeratorAbs, periodAbs) = (AbsI(numerator), AbsI(period));

        return (periodAbs * currentDivisor) / GreatestCommonDivisorI(currentDivisor, periodAbs);
    }

    function MaybeFactorsFromPeriod(modulus : Int, generator : Int, period : Int) 
    : (Bool, (Int, Int)) {
        if (period % 2 == 0) {

            let halfPower = ExpModI(generator, period / 2, modulus);

            if (halfPower != modulus - 1) {

                let factor = MaxI(
                    GreatestCommonDivisorI(halfPower - 1, modulus), 
                    GreatestCommonDivisorI(halfPower + 1, modulus)
                );
                
                return (true, (factor, modulus / factor));
            } else {
                return (false, (1,1));
            }
        } else {

            Message("Estimated period was odd, trying again.");
            return (false, (1,1));
        }
    }
}
