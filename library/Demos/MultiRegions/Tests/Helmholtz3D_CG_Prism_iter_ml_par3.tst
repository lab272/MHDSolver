<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 3D CG, prisms, Neumann BCs, iterative ML, Par(3)</description>
    <executable>Helmholtz3D</executable>
    <parameters>--use-scotch -I WeightPartitions=Boundary Helmholtz3D_Prism.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Helmholtz3D_Prism.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="3e-5">0.000542783</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">0.0102924</value>
        </metric>
    </metrics>
</test>
