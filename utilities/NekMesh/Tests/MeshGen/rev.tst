<NEKTAR>
    <MESHING>

        <INFORMATION>
            <I PROPERTY="CADFile" VALUE="sphere.STEP" />
            <I PROPERTY="MeshType" VALUE="3D" />
        </INFORMATION>

        <PARAMETERS>
            <P PARAM="MinDelta" VALUE="0.01"    />
            <P PARAM="MaxDelta" VALUE="0.2"     />
            <P PARAM="EPS"      VALUE="0.01"     />

            <P PARAM="Order"    VALUE="4"       />
        </PARAMETERS>

        <BOOLPARAMETERS>
            <P VALUE="SurfaceOptimiser" />
        </BOOLPARAMETERS>

    </MESHING>
</NEKTAR>
