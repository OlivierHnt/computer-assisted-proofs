# Problems declared and not written

struct ManifoldProblem{Ff,FDf,P,C,L,V,BX} <: AbstractSequenceProblem
    f      :: Ff
    Df     :: FDf
    p      :: P
    c      :: C
    lambda :: L
    V      :: V
    X      :: BX
end

struct ConnectingOrbitProblem{Ff,FDf,P,U,S,T,BX} <: AbstractSequenceProblem
    f        :: Ff
    Df       :: FDf
    p        :: P
    unstable :: U
    stable   :: S
    tau      :: T
    X        :: BX
end
