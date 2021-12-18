namespace Microsoft.Quantum.Sample.Teleportation {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;

  
    operation Teleport (msg : Qubit, target : Qubit) : Unit {
        use register = Qubit();
        H(register);
        CNOT(register, target);

        CNOT(msg, register);
        H(msg);

        if (MResetZ(msg) == One) { Z(target); }
        if (IsResultOne(MResetZ(register))) { X(target); }
    }

    operation TeleportClassicalMessage (message : Bool) : Bool {
        use (msg, target) = (Qubit(), Qubit());

        if (message) {
            X(msg);
        }

        Teleport(msg, target);

        return MResetZ(target) == One;
    }

    operation TeleportRandomMessage () : Unit {
        use (msg, target) = (Qubit(), Qubit());
        PrepareRandomMessage(msg);

        Teleport(msg, target);

        if (MeasureIsPlus(target))  { Message("Received |+>"); }
        if (MeasureIsMinus(target)) { Message("Received |->"); }

        Reset(msg);
        Reset(target);
    }
}



