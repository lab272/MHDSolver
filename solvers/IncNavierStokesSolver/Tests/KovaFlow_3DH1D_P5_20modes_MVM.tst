<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes (MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_MVM.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.95163e-05</value>
            <value variable="v" tolerance="1e-12">7.6401e-06</value>
            <value variable="w" tolerance="1e-12">1.33702e-05</value>
	    <value variable="p" tolerance="1e-12">6.04171e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000106884</value>
            <value variable="v" tolerance="1e-12">3.9563e-05</value>
            <value variable="w" tolerance="1e-12">3.9394e-05</value>
	    <value variable="p" tolerance="1e-12">0.000513613</value>
        </metric>
    </metrics>
</test>