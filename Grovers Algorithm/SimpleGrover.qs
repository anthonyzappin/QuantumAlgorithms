namespace Microsoft.Quantum.Samples.SimpleGrover {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;

    @EntryPoint()
    operation SearchForMarkedInput(nQubits : Int) : Result[] {
        use qubits = Qubit[nQubits];
        PrepareUniform(qubits);
        for idxIteration in 0..NIterations(nQubits) - 1 {
            ReflectAboutMarked(qubits);
            ReflectAboutUniform(qubits);
        }
        return ForEach(MResetZ, qubits);
    }

    function NIterations(nQubits : Int) : Int {
        let nItems = 1 <<< nQubits; 
        let angle = ArcSin(1. / Sqrt(IntAsDouble(nItems)));
        let nIterations = Round(0.25 * PI() / angle - 0.5);
        return nIterations;
    }

}
