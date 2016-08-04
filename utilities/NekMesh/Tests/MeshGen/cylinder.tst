<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Simple geometry with cylinder</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list cylinder.mcf cylinder.xml:xml:test</parameters>
    <files>
        <file description="Input File">cylinder.mcf</file>
        <file description="Input File 2">cylinder.STEP</file>
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
