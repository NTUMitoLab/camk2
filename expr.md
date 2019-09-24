# ODE functions
```math
\begin{aligned}
\frac{dCaM}{dt} =&  - CaM \cdot k1p \cdot ca + k1m \cdot CaMCa1 + \frac{kdissoCa \cdot kmCaM^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKII_{CaMCa4} + \frac{kdissoCa2 \cdot kmCaM^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKIIp_{CaMCa4} \\
\frac{dCaMCa1}{dt} =& CaM \cdot k1p \cdot ca - k1m \cdot CaMCa1 - CaMCa1 \cdot k2p \cdot ca + k2m \cdot CaMCa2 \\
\frac{dCaMCa2}{dt} =& CaMCa1 \cdot k2p \cdot ca - k2m \cdot CaMCa2 - CaMCa2 \cdot k3p \cdot ca + k3m \cdot CaMCa3 \\
\frac{dCaMCa3}{dt} =& CaMCa2 \cdot k3p \cdot ca - k3m \cdot CaMCa3 - CaMCa3 \cdot k4p \cdot ca + k4m \cdot CaMCa4 \\
\frac{dCaMCa4}{dt} =& CaMCa3 \cdot k4p \cdot ca - k4m \cdot CaMCa4 - kasso \cdot CaMCa4 \cdot CaMKII + \frac{kdisso \cdot ca^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKII_{CaMCa4} + \frac{kdisso2 \cdot ca^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKIIp_{CaMCa4} \\
\frac{dCaMKII}{dt} =&  - kasso \cdot CaMCa4 \cdot CaMKII + \frac{kdisso \cdot ca^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKII_{CaMCa4} + \frac{kdissoCa \cdot kmCaM^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKII_{CaMCa4} + \frac{kcatPP1 \cdot CaMKIIp}{KmPP1 + CaMKIIp} \cdot PP1 \\
\frac{dCaMKII_{CaMCa4}}{dt} =& kasso \cdot CaMCa4 \cdot CaMKII - \frac{kdisso \cdot ca^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKII_{CaMCa4} - \frac{kdissoCa \cdot kmCaM^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKII_{CaMCa4} - CaMKII_{CaMCa4} \cdot \frac{kcat \cdot ATP}{KmATP + ATP} \cdot \left( 1 - \left( \frac{CaMKII}{\Sigma CaMKII} \right)^{2} \right) + \frac{kcatPP1 \cdot CaMKIIp_{CaMCa4}}{KmPP1 + CaMKIIp_{CaMCa4}} \cdot PP1 \\
\frac{dCaMKIIp_{CaMCa4}}{dt} =&  - \frac{kdisso2 \cdot ca^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKIIp_{CaMCa4} - \frac{kdissoCa2 \cdot kmCaM^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKIIp_{CaMCa4} + CaMKII_{CaMCa4} \cdot \frac{kcat \cdot ATP}{KmATP + ATP} \cdot \left( 1 - \left( \frac{CaMKII}{\Sigma CaMKII} \right)^{2} \right) - \frac{kcatPP1 \cdot CaMKIIp_{CaMCa4}}{KmPP1 + CaMKIIp_{CaMCa4}} \cdot PP1 \\
\frac{dCaMKIIp}{dt} =& \frac{kdisso2 \cdot ca^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKIIp_{CaMCa4} + \frac{kdissoCa2 \cdot kmCaM^{3}}{kmCaM^{3} + ca^{3}} \cdot CaMKIIp_{CaMCa4} - \frac{kcatPP1 \cdot CaMKIIp}{KmPP1 + CaMKIIp} \cdot PP1
\end{aligned}
```
