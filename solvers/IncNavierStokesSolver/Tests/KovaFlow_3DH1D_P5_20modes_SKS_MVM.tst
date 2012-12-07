<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes - Skew-Symmetric advection(MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.00306e-05</value>
            <value variable="v" tolerance="1e-12">7.82855e-06</value>
            <value variable="w" tolerance="1e-12">1.35226e-05</value>
	    <value variable="p" tolerance="1e-12">6.29623e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000109399</value>
            <value variable="v" tolerance="1e-12">3.96664e-05</value>
            <value variable="w" tolerance="1e-12">4.0213e-05</value>
	    <value variable="p" tolerance="1e-12">0.000565171</value>
        </metric>
    </metrics>
</test>