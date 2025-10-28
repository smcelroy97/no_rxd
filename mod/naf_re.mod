TITLE naf_re.mod
 
COMMENT
 Fast sodium current for RE and TC cells that was supposed to be defined in appendix  Bazhenov 1998 (doi: https://doi.org/10.1152/jn.1998.79.5.2730), and referenced in 
 Bazhenov et. al. 2002 (J Neuro) and Chen et. al. 2012 (J. Physiol.)
 In reality, this code apes that from the C++ code accompanying Bazhenov 2002
 (found here: https://modeldb.science/28189?tab=1) 
 See the code for class 'INaK', 
 and then also--crucially--the change of parameters around line 1500 of Bazhenov 2002 neur271.cell
 (especially for Vtr)
 This code is adapted from hh.mod distributed with NEURON source code
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX naf_re
        USEION na READ ena WRITE ina
        RANGE gnabar, gna
        RANGE minf, hinf, mtau, htau
	:THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = 0.0015 (S/cm2)	<0,1e9>
        ena = 50 (mV)
}
 
STATE {
        m h
}
 
ASSIGNED {
        v (mV)
        :celsius (degC)

		gna (S/cm2)
        ina (mA/cm2)
        minf hinf
		mtau (ms) htau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
		ina = gna*(v - ena) 
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
}
 
? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  a, b, sum
UNITSOFF
		:"m" sodium activation system
        a = 0.32 * vtrap(-(v+37),4) :for TC cell, it should be 27 (instead of 37) in order to match C++ code (around line 1500 of neur271.c)
        b =  0.28 * vtrap(v+10,5) :for TC cell, it should be 0 (instead of 10) in order to match C++ code (around line 1500 of neur271.c)
        sum = a + b
		mtau = 1.0/sum
        minf = a/sum
                :"h" sodium inactivation system
        
        a = 0.128*exp(-(v+33)/18) :for TC cell, it should be 23 (instead of 33) in order to match C++ code (around line 1500 of neur271.c)
        b = 4/(exp(-(v+10)/5) + 1) :for TC cell, it should be 0 (instead of 10) in order to match C++ code (around line 1500 of neur271.c)
        sum = a + b
		htau = 1.0/sum
		hinf = a/sum
        
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2) :derived by Taylor expanding exp to quadratic term, then binomial approximation
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON


