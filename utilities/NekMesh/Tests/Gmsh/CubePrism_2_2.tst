<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v2.2) high-order prism cube, order 2</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list CubePrism_2_2.msh CubePrism_2_2.xml:xml:test</parameters>
    <files>
        <file description="Input File">CubePrism_2_2.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
