#Dr. Smay helped me build this code, he also provided the starting template
# Ashley helped when making this code
# ChatGPT helped with this code

#region imports
import sys
from ThermoStateCalc_update import Ui__frm_StateCalculator
from pyXSteam.XSteam import XSteam
from PyQt5.QtWidgets import QWidget, QApplication
from UnitConversion import UC
from scipy.optimize import fsolve
#endregion

#region class definitions
class state:
    def __init__(self):
        """
        Initializes a new instance of the state class, representing a thermodynamic state.

        Attributes:
            t (float, optional): Temperature of the thermodynamic state, initially None.
            v (float, optional): Specific volume of the thermodynamic state, initially None.
            p (float, optional): Pressure of the thermodynamic state, initially None.
            h (float, optional): Enthalpy of the thermodynamic state, initially None.
            u (float, optional): Internal energy of the thermodynamic state, initially None.
            s (float, optional): Entropy of the thermodynamic state, initially None.
            x (float, optional): Quality (vapor fraction) of the thermodynamic state, initially None.
                                 This is relevant in two-phase regions where 0 <= x <= 1.
        """
        self.t= None
        self.v= None
        self.p= None
        self.h= None
        self.u= None
        self.s= None
        self.x= None
class state2:
    def __int__(self):
        self.t = None
        self.v = None
        self.p = None
        self.h = None
        self.u = None
        self.s = None
        self.x = None
class main_window(QWidget,Ui__frm_StateCalculator):
    def __init__(self):
        """
        Initializes the main window of the application, setting up the user interface,
        connecting signals and slots, and initializing the steam table with the default
        unit system.

        The constructor performs the following actions:
        - Calls the constructor of the superclass to initialize the inherited QWidget.
        - Sets up the user interface as defined in Ui__frm_StateCalculator.
        - Connects the widget signals to their respective slots to handle user interactions.
        - Initializes the steam table object with the Metric (MKS) unit system.
        - Sets the initial units for the application to SI units.
        - Calls the method to set the units for the displayed properties.
        - Shows the main window on the screen
        """
        super().__init__()
        self.setupUi(self)
        self.SetupSlotsAndSignals()
        self.steamTable=XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.currentUnits='SI'
        self.setUnits()
        self.show()

    def SetupSlotsAndSignals(self):
        """
        Connects the UI elements' signals to the respective slots to handle user interactions.

    This method establishes the connections between the interactive UI components
    (like radio buttons, combo boxes, and buttons) and their corresponding methods
    to ensure the application responds to user actions appropriately. Specifically, it:

    - Connects the radio buttons for unit selection (_rdo_English, _rdo_SI) to the setUnits method,
      which updates the units displayed in the application.
    - Connects the combo boxes for property selection in both states (_cmb_Property1_1, _cmb_Property2_1
      for state 1 and _cmb_Property1_2, _cmb_Property2_2 for state 2) to the setUnits method, which
      will update the units based on the selected properties.
    - Connects the calculate button (_pb_Calculate) to the onCalculate method, which triggers the
      calculation of the thermodynamic properties based on the user's input.
        :return:None
        """

        self._rdo_English.clicked.connect(self.setUnits)
        self._rdo_SI.clicked.connect(self.setUnits)
        self._cmb_Property1_1.currentIndexChanged.connect(self.setUnits) # This is state 1
        self._cmb_Property2_1.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property1_2.currentIndexChanged.connect(self.setUnits) # This is state 2
        self._cmb_Property2_2.currentIndexChanged.connect(self.setUnits)
        self._pb_Calculate.clicked.connect(self.onCalculate)
        pass

    def onCalculate(self):
        """
        Handles the calculation process when the calculate button is clicked.

    This method orchestrates the calculation of thermodynamic properties for two states
    and the change between these states. It performs the following actions:

    - Calls `calculateProperties` to compute the properties of the first state based on
      the user inputs from the UI and assigns the result to `self.state1`.
    - Calls `calculateProperties2` to compute the properties of the second state based
      on the user inputs from the UI and assigns the result to `self.state2`.
    - Calls `calculateDelta` with `self.state1` and `self.state2` to calculate the change
      in properties between the two states and assigns the result to `self.stateChange`.

    These steps collectively update the application's state to reflect the new calculations
    based on the user's input, showing the individual state properties and their differences.
        :return:
        """
        self.state1 = self.calculateProperties()
        self.state2 = self.calculateProperties2()
        self.stateChange = self.calculateDelta(self.state1, self.state2)

    def setUnits(self):
        """
        This sets the units for the selected specified properties.
        Units for the thermodynamic properties are set upon pushing calculate button.
        I added part of the same code again to account for state 2.
        :return: None
        """
        #set the units system based on selected radio button
        #also, determine if a units change is required
        SI=self._rdo_SI.isChecked()
        newUnits='SI' if SI else 'EN'
        UnitChange = self.currentUnits != newUnits  # compare new units to current units
        self.currentUnits = newUnits

        if SI:
            self.steamTable=XSteam(XSteam.UNIT_SYSTEM_MKS)
            self.l_Units = "m"
            self.p_Units = "bar"
            self.t_Units = "C"
            self.m_Units = "kg"
            self.time_Units = "s"
            self.energy_Units = "W"
            self.u_Units = "kJ/kg"
            self.h_Units = "kJ/kg"
            self.s_Units = "kJ/kg*C"
            self.v_Units = "m^3/kg"
        else:
            self.steamTable=XSteam(XSteam.UNIT_SYSTEM_FLS)
            self.l_Units = "ft"
            self.p_Units = "psi"
            self.t_Units = "F"
            self.m_Units = "lb"
            self.time_Units = "s"
            self.energy_Units = "btu"
            self.u_Units = "btu/lb"
            self.h_Units = "btu/lb"
            self.s_Units = "btu/lb*F"
            self.v_Units = "ft^3/lb"

        #read selected Specified Properties from combo boxes
        SpecifiedProperty1_1 = self._cmb_Property1_1.currentText() #Property 1 state 1
        SpecifiedProperty2_1 = self._cmb_Property2_1.currentText() #Property 2 state 1
        SpecifiedProperty1_2 = self._cmb_Property1_2.currentText() #Property 1 state 2
        SpecifiedProperty2_2 = self._cmb_Property2_2.currentText() #property 2 state 2
        #read numerical values for selected properties
        SP=[float(self._le_Property1_1.text()), float(self._le_Property2_1.text()), float(self._le_Property2_2.text()),float(self._le_Property2_2.text())]

        #set units labels and convert values if needed
        if 'Pressure' in SpecifiedProperty1_1:
            self._lbl_Property1_1_Units.setText(self.p_Units)
            if UnitChange:  # note that I only should convert if needed.  Not if I double click on SI or English
                SP[0]=SP[0]*UC.psi_to_bar if SI else SP[0]*UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty1_1:
            self._lbl_Property1_1_Units.setText(self.t_Units)
            if UnitChange:
                SP[0] = UC.F_to_C(SP[0]) if SI else UC.C_to_F(SP[0])
        elif 'Energy' in SpecifiedProperty1_1:
            self._lbl_Property1_1_Units.setText(self.u_Units)
            if UnitChange:
                SP[0]=SP[0]*UC.btuperlb_to_kJperkg if SI else SP[0]*UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty1_1:
            self._lbl_Property1_1_Units.setText(self.h_Units)
            if UnitChange:
                SP[0] = SP[0] * UC.btuperlb_to_kJperkg if SI else SP[0] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty1_1:
            self._lbl_Property1_1_Units.setText(self.s_Units)
            if UnitChange:
                SP[0] = SP[0] * UC.btuperlbF_to_kJperkgC if SI else SP[0] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty1_1:
            self._lbl_Property1_1_Units.setText(self.v_Units)
            if UnitChange:
                SP[0]=SP[0]*UC.ft3perlb_to_m3perkg if SI else SP[0]*UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty1_1:
            self._lbl_Property1_1_Units.setText("")

        if 'Pressure' in SpecifiedProperty2_1:
            self._lbl_Property2_1_Units.setText(self.p_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.psi_to_bar if SI else SP[1] * UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty2_1:
            self._lbl_Property2_1_Units.setText(self.t_Units)
            if UnitChange:
                SP[1] = UC.F_to_C(SP[1]) if SI else UC.C_to_F(SP[1])
        elif 'Energy' in SpecifiedProperty2_1:
            self._lbl_Property2_1_Units.setText(self.u_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.btuperlb_to_kJperkg if SI else SP[1] * UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty2_1:
            self._lbl_Property2_1_Units.setText(self.h_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.btuperlb_to_kJperkg if SI else SP[1] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty2_1:
            self._lbl_Property2_1_Units.setText(self.s_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.btuperlbF_to_kJperkgC if SI else SP[1] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty2_1:
            self._lbl_Property2_1_Units.setText(self.v_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.ft3perlb_to_m3perkg if SI else SP[1] * UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty2_1:
            self._lbl_Property2_1_Units.setText("")

        if 'Pressure' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.p_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.psi_to_bar if SI else SP[2] * UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.t_Units)
            if UnitChange:
                SP[2] = UC.F_to_C(SP[2]) if SI else UC.C_to_F(SP[2])
        elif 'Energy' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.u_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.btuperlb_to_kJperkg if SI else SP[2] * UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.h_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.btuperlb_to_kJperkg if SI else SP[2] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.s_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.btuperlbF_to_kJperkgC if SI else SP[2] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.v_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.ft3perlb_to_m3perkg if SI else SP[2] * UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText("")

        if 'Pressure' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.p_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.psi_to_bar if SI else SP[3] * UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.t_Units)
            if UnitChange:
                SP[3] = UC.F_to_C(SP[3]) if SI else UC.C_to_F(SP[3])
        elif 'Energy' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.u_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.btuperlb_to_kJperkg if SI else SP[3] * UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.h_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.btuperlb_to_kJperkg if SI else SP[3] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.s_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.btuperlbF_to_kJperkgC if SI else SP[3] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.v_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.ft3perlb_to_m3perkg if SI else SP[3] * UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText("")

        self._le_Property1_1.setText("{:0.3f}".format(SP[0]))
        self._le_Property2_1.setText("{:0.3f}".format(SP[1]))
        self._le_Property1_2.setText("{:0.3f}".format(SP[2]))
        self._le_Property2_2.setText("{:0.3f}".format(SP[3]))



    def clamp(self, x, low, high):
        """
        This clamps a float x between a high and low limit inclusive
        :param x:The numeric value to be clamped
        :param low: The lower bound of the range
        :param high: The higher bound of the range
        :return: x
        """
        if x<low:
            return low
        if x>high:
            return high
        return x

    def between(self, x, low, high):
        """
        Tells if x is between low and high inclusive
        :param x: The numeric value to be checked
        :param low: The lower bound of the range
        :param high: The upper bound of the range
        :return: Bool: True if 'x' is between 'low' and 'high' otherwise False.
        """
        if x>=low and x<=high:
            return True
        return False

    def getSatProps_p(self, p):
        """
        Given a pressure, calculate the saturated properties for that isobar
        :param p: The pressure at which to calculate saturated properties
        :return: None
        """
        self.tSat=self.steamTable.tsat_p(p)
        self.pSat=p
        self.vf=self.steamTable.vL_p(p)
        self.vg=self.steamTable.vV_p(p)
        self.hf=self.steamTable.hL_p(p)
        self.hg=self.steamTable.hV_p(p)
        self.uf=self.steamTable.uL_p(p)
        self.ug=self.steamTable.uV_p(p)
        self.sf=self.steamTable.sL_p(p)
        self.sg=self.steamTable.sV_p(p)
        self.vgf=self.vg-self.vf
        self.hgf=self.hg-self.hf
        self.sgf=self.sg-self.sf
        self.ugf=self.ug-self.uf

    def getSatProps_t(self, t):
        """
        Given a temperature, calculate the saturation pressure and then
        calculate all other saturated properties
        :param t: The temperature at which to calculate the saturated properties, typically in degrees Celsius or Kelvin
        :return: None
        """
        self.tSat=t
        self.pSat=self.steamTable.psat_t(t)
        self.getSatProps_p(self.pSat)

    def makeLabel_1Phase(self):
        """
        Given that I have T&P, find the other properties and make the label. This is for state 1
        :return:None
        """
        self._lbl_State_1.setText("Region = {:}".format(self.region))
        stProps = "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
        stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat,
                                                                               self.t_Units)
        self.u = self.steamTable.u_pt(self.p, self.t)
        self.h = self.steamTable.h_pt(self.p, self.t)
        self.s = self.steamTable.s_pt(self.p, self.t)
        self.v = self.steamTable.v_pt(self.p, self.t)
        self.x = 1.0 if self.t > self.steamTable.tsat_p(self.p) else 0.0
        stProps += "\nInternal Energy = {:0.3f} ({:})".format(self.u, self.u_Units)
        stProps += "\nEnthalpy = {:0.3f} ({:})".format(self.h, self.h_Units)
        stProps += "\nEntropy = {:0.3f} ({:})".format(self.s, self.s_Units)
        stProps += "\nSpecific Volume = {:0.3f} ({:})".format(self.v, self.v_Units)
        stProps += "\nQuality = {:0.3f}".format(self.x)
        self.stProps=stProps

    def makeLabel_1Phase_state2(self):
        """
        Given that I have T&P, find the other properties and make the label.This is for state 2
        :return: None
        """
        self._lbl_State_2.setText("Region = {:}".format(self.region))
        stProps = "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
        stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat,
                                                                               self.t_Units)
        self.u = self.steamTable.u_pt(self.p, self.t)
        self.h = self.steamTable.h_pt(self.p, self.t)
        self.s = self.steamTable.s_pt(self.p, self.t)
        self.v = self.steamTable.v_pt(self.p, self.t)
        self.x = 1.0 if self.t > self.steamTable.tsat_p(self.p) else 0.0
        stProps += "\nInternal Energy = {:0.3f} ({:})".format(self.u, self.u_Units)
        stProps += "\nEnthalpy = {:0.3f} ({:})".format(self.h, self.h_Units)
        stProps += "\nEntropy = {:0.3f} ({:})".format(self.s, self.s_Units)
        stProps += "\nSpecific Volume = {:0.3f} ({:})".format(self.v, self.v_Units)
        stProps += "\nQuality = {:0.3f}".format(self.x)
        self.stProps=stProps

    def makeLabel_2Phase(self):
        """
        Given P and x, find all other properties and make the label. This is for state 1
        :return: None
        """
        self._lbl_State_1.setText("Region = {:}".format(self.region))
        stProps = "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
        stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat,
                                                                               self.t_Units)
        stProps += "\nInternal Energy = u={:0.3f} ({:})".format(self.uf + self.x * self.ugf, self.u_Units)
        stProps += "\nEnthalpy = h={:0.3f} ({:})".format(self.hf + self.x * self.hgf, self.h_Units)
        stProps += "\nEntropy = s={:0.3f} ({:})".format(self.sf + self.x * self.sgf, self.s_Units)
        stProps += "\nSpecific Volume = v={:0.5f} ({:})".format(self.vf + self.x * self.vgf, self.v_Units)
        stProps += "\nQuality = {:0.3f}".format(self.x)
        self.stProps=stProps

    def makeLabel_2Phase_state2(self):
        """
        Given P and x, find all other properties and make the label. This is for state 2
        :return: None
        """
        self._lbl_State_2.setText("Region = {:}".format(self.region))
        stProps = "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
        stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat,
                                                                               self.t_Units)
        stProps += "\nInternal Energy = u={:0.3f} ({:})".format(self.uf + self.x * self.ugf, self.u_Units)
        stProps += "\nEnthalpy = h={:0.3f} ({:})".format(self.hf + self.x * self.hgf, self.h_Units)
        stProps += "\nEntropy = s={:0.3f} ({:})".format(self.sf + self.x * self.sgf, self.s_Units)
        stProps += "\nSpecific Volume = v={:0.5f} ({:})".format(self.vf + self.x * self.vgf, self.v_Units)
        stProps += "\nQuality = {:0.3f}".format(self.x)
        self.stProps=stProps

    def calculateProperties(self):
        """
        Calculates the thermodynamic state variables based on specified values.
        I have thermodynamic variables:  P, T, v, h, u, s and x (7 things) from which I am choosing two.
        Possible number of permutations:  7!/5! =42.
        But, order of the two things does not matter, so 42/2=21
        PT, Pv, Ph, Pu, Ps, Px (6)
        Tv, Th, Tu, Ts, Tx (5)
        vh, vu, vs, vx (4)
        hu, hs, hx (3)
        us, ux (2)
        sx (1)
        Total of 21 cases to deal with.  I will attack them in the order shown above
        This part is only calculating the properties for state 1. It was copied and used for Calculate
        properties2 to calculate the properties for state 2.
        :return: nothing
        """
        # Step 1: read which properties are being specified from the combo boxes
        SP=[self._cmb_Property1_1.currentText()[-2:-1], self._cmb_Property2_1.currentText()[-2:-1]]
        if SP[0]==SP[1]:
            self._lbl_Warning.setText("Warning:  You cannot specify the same property twice.")
        else:
            self._lbl_Warning.setText("")
        SP[0]=SP[0].lower()
        SP[1]=SP[1].lower()

        # Step 2: select the proper case from the 21.  Note that PT is the same as TP etc.
        if SP[0]=='p' or SP[1]=='p':
            oFlipped= SP[0] != 'p'
            SP1 = SP[0] if oFlipped else SP[1]
        #case 1:  pt or tp
            if SP1=='t':
                f1=float(self._le_Property1_1.text())
                f2=float(self._le_Property2_1.text())
                self.p= f1 if not oFlipped else f2
                self.t= f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.tSat=round(self.tSat,3)  # I will compare at 3 three decimal places
                # compare T to TSat
                if self.t<self.tSat or self.t>self.tSat:
                    self.region = "sub-cooled liquid" if self.t<self.tSat else "super-heated vapor"
                    self.makeLabel_1Phase()
                else:  #this is ambiguous since at saturated temperature
                    self.region = "two-phase"
                    self.stProps ="Region = {:}".format(self.region)
                    self.stProps += "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
                    self.stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat, self.t_Units)
                    self.stProps += "\nInternal Energy = uf={:0.3f}, ug={:0.3f} ({:})".format(self.uf, self.ug, self.u_Units)
                    self.stProps += "\nEnthalpy = hf={:0.3f}, hg={:0.3f} ({:})".format(self.hf,self.hg, self.h_Units)
                    self.stProps += "\nEntropy = sf={:0.3f}, sg={:0.3f} ({:})".format(self.sf, self.sg, self.s_Units)
                    self.stProps += "\nSpecific Volume = vf={:0.3f}, vg={:0.3f} ({:})".format(self.vf, self.vg, self.v_Units)
                    self.stProps += "\nQuality = unknown"
        #case 2: pv or vp
            elif SP1=='v':
                f1=float(self._le_Property1_1.text())
                f2=float(self._le_Property2_1.text())
                self.p= f1 if not oFlipped else f2
                self.v= f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.vf=round(self.vf,5)
                self.vg=round(self.vg,3)
                # compare v to vf and vg
                if self.v<self.vf or self.v>self.vg:
                    self.region = "sub-cooled liquid" if self.v<self.vf else "super-heated vapor"
                    #since I can't find properties using v, I will use fsolve to find T
                    dt=1.0 if self.v>self.vg else -1.0
                    fn1 = lambda T: self.v - self.steamTable.v_pt(self.p, T)
                    self.t = fsolve(fn1, [self.tSat + dt])[0]
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.v-self.vf)/(self.vgf)
                    self.t=self.tSat
                    self.makeLabel_2Phase()
        #case 3 pu or up
            elif SP1=='u':
                f1=float(self._le_Property1_1.text())
                f2=float(self._le_Property2_1.text())
                self.p= f1 if not oFlipped else f2;
                self.u= f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u<self.uf or self.u>self.ug:
                    self.region = "sub-cooled liquid" if self.u<self.uf else "super-heated vapor"
                    #since I can't find properties using u, I will use fsolve to find T
                    dt=1.0 if self.u>self.ug else -1.0
                    fn3 = lambda T: self.u - self.steamTable.u_pt(self.p, T)
                    self.t = fsolve(fn3, [self.tSat + dt])[0]
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.u-self.uf)/(self.ugf)
                    self.t=self.tSat
                    self.makeLabel_2Phase()
        #case 4 ph or hp
            elif SP1 == 'h':
                f1=float(self._le_Property1_1.text())
                f2=float(self._le_Property2_1.text())
                self.p= f1 if not oFlipped else f2;
                self.h= f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare h to hf and hg
                if self.h < self.hf or self.h > self.hg:
                    self.region = "sub-cooled liquid" if self.h < self.hf else "super-heated vapor"
                    self.t=self.steamTable.t_ph(self.p, self.h)
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.h - self.hf) / (self.hgf)
                    self.t=self.tSat
                    self.makeLabel_2Phase()
        #case 5 ps or sp
            elif SP1 == 's':
                f1=float(self._le_Property1_1.text())
                f2=float(self._le_Property2_1.text())
                self.p= f1 if not oFlipped else f2;
                self.s= f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    self.t=self.steamTable.t_ps(self.p, self.s)
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
        #case 6 px or xp
            elif SP1 == 'x':
                self.region="two-phase"
                f1=float(self._le_Property1_1.text())
                f2=float(self._le_Property2_1.text())
                self.p= f1 if not oFlipped else f2
                self.x= f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.t=self.tSat
                self.x=self.clamp(self.x,0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0]=='t' or SP[1]=='t':
            oFlipped = SP[0] != 't'
            SP1 = SP[0] if oFlipped else SP[1]
        #case 7:  tv or vt
            if SP1=='v':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.t = f1 if not oFlipped else f2
                self.v = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                self.vf=round(self.vf,5)
                self.vg=round(self.vg,3)
                # compare v to vf and vg
                if self.v<self.vf or self.v>self.vg:
                    self.region = "sub-cooled liquid" if self.v<self.vf else "super-heated vapor"
                    #since I can't find properties using v, I will use fsolve to find P
                    dp=-0.1 if self.v>self.vg else 0.1
                    fn3 = lambda P: [self.v - self.steamTable.v_pt(P, self.t)]
                    self.p = fsolve(fn3, [self.pSat + dp])[0]
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.v-self.vf)/(self.vgf)
                    self.p=self.pSat
                    self.makeLabel_2Phase()
        #case 8:  tu or ut
            elif SP1=='u':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.t = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u<self.uf else "super-heated vapor"
                    #since I can't find properties using u, I will use fsolve to find P
                    dp=0.1 if self.u>self.ug else -0.1
                    fn8 = lambda P: self.u - self.steamTable.u_pt(self.t, P)
                    self.p = fsolve(fn8, [self.pSat + dp])[0]
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.u-self.uf)/(self.ugf)
                    self.p=self.pSat
                    self.makeLabel_2Phase()
        #case 9:  th or ht
            elif SP1 == 'h':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.t = f1 if not oFlipped else f2
                self.h = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare h to hf and hg
                if self.h < self.hf or self.h > self.hg:
                    self.region = "sub-cooled liquid" if self.h < self.hf else "super-heated vapor"
                    self.p=self.steamTable.p_th(self.t, self.h)
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p=self.pSat
                    self.x = (self.h - self.hf) / (self.hgf)
                    self.makeLabel_2Phase()
        #case 10:  ts or st
            elif SP1 == 's':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.t = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    self.p=self.steamTable.p_ts(self.t, self.s)
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p=self.pSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
        #case 11:  tx or xt
            elif SP1 == 'x':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.t = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region="two-phase"
                self.getSatProps_t(self.t)
                self.p=self.pSat
                self.x = float(self._le_Property2_1.text())
                self.x=self.clamp(self.x,0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0] == 'v' or SP[1] == 'v':
            oFlipped = SP[0] != 'v'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 12:  vh or hv
            if SP1 == 'h':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.h = f2 if not oFlipped else f1
                def fn12(P):
                    # could be single phase or two-phase, but both v&h have to match at same x
                    self.getSatProps_p(P)
                    if self.between(self.h,self.hf, self.hg):
                        self.x=(self.h-self.hf)/self.hgf
                        return self.v-(self.vf+self.x*self.vgf)
                    # could be single phase
                    return self.v-self.steamTable.v_ph(P,self.h)
                self.p=fsolve(fn12,[1.0])[0]
                self.t=self.steamTable.t_ph(self.p, self.h)
                self.getSatProps_p(self.p)
                self.vf = round(self.vf, 5)
                self.vg = round(self.vg, 3)
                # compare v to vf and vg
                if self.v < self.vf or self.v > self.vg:
                    self.region = "sub-cooled liquid" if self.v < self.vf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.v - self.vf) / (self.vgf)
                    self.makeLabel_2Phase()
            # case 13:  vu or uv
            elif SP1 == 'u':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.v = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                # use fsolve to fing P&T at this v & u
                def fn13(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.u,self.uf, self.ug):
                        self.t=self.tSat
                        self.x=(self.u-self.uf)/self.ugf
                        return [self.v-(self.vf+self.x*self.vgf),0]
                    return [self.v-self.steamTable.v_pt(PT[0],PT[1]),self.u-self.steamTable.u_pt(PT[0],PT[1])]
                props=fsolve(fn13, [1,100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.p = self.pSat
                    self.t=self.tSat
                    self.makeLabel_2Phase()
            # case 14:  vs os sv
            elif SP1 == 's':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.v = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                def fn14(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x=(self.s-self.sf)/self.sgf
                        return [self.v-self.vf-self.x*self.vgf, 0.0]
                    return [self.v - self.steamTable.v_pt(PT[0], PT[1]),
                                  self.s - self.steamTable.s_pt(PT[0], PT[1])]
                props = fsolve(fn14, [1, 100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
            # case 15:  vx or xv
            elif SP1 == 'x':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.v = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn15(p):
                    self.getSatProps_p(p)
                    return self.v -(self.vf+self.x*self.vgf)
                props=fsolve(fn15, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_1.text()),0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0] == 'h' or SP[1] == 'h':
            oFlipped = SP[0] != 'h'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 16:  hu or uh
            if SP1 == 'u':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.h = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                # use fsolve to fing P&T at this v & u
                def fn16(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.u, self.uf, self.ug):
                        self.x=(self.u-self.uf)/self.ugf
                        return [self.h-self.hf-self.x*self.hgf, 0.0]
                    return [self.h-self.steamTable.h_pt(PT[0],PT[1]),self.u-self.steamTable.u_pt(PT[0],PT[1])]
                props=fsolve(fn16, [1,100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.p = self.pSat
                    self.t=self.tSat
                    self.makeLabel_2Phase()
            # case 17:  hs or sh
            elif SP1 == 's':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.h = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                def fn17(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x=(self.s-self.sf)/self.sgf
                        return [self.h-self.hf-self.x*self.hgf, 0.0]
                    return [self.h-self.steamTable.h_pt(PT[0],PT[1]),self.s-self.steamTable.s_pt(PT[0],PT[1])]
                props = fsolve(fn17, [1, 100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
            # case 18:  hx or xh
            elif SP1 == 'x':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.v = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn18(p):
                    self.getSatProps_p(p)
                    return self.h -(self.hf+self.x*self.hgf)
                props=fsolve(fn18, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_1.text()),0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0] == 'u' or SP[1] == 'u':
            oFlipped = SP[0] != 'u'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 19:  us or su
            if SP1 == 's':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.u = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                def fn19(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x = (self.s - self.sf) / self.sgf
                        return [self.u - self.uf - self.x * self.ugf, 0.0]
                    return [self.u - self.steamTable.u_pt(PT[0], PT[1]),
                            self.s - self.steamTable.s_pt(PT[0], PT[1])]
                props = fsolve(fn19, [1, 100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t=self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
            # case 20:  ux or xu
            elif SP1 == 'x':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.u = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn20(p):
                    self.getSatProps_p(p)
                    return self.h -(self.hf+self.x*self.hgf)
                props=fsolve(fn20, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_1.text()),0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0] == 's' or SP[1] == 's':
            oFlipped = SP[0] != 's'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 21:  sx or xs
            if SP1 == 'x':
                f1 = float(self._le_Property1_1.text())
                f2 = float(self._le_Property2_1.text())
                self.s = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn21(p):
                    self.getSatProps_p(p)
                    return self.h -(self.hf+self.x*self.hgf)
                props=fsolve(fn21, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_1.text()),0.0,1.0)
                self.makeLabel_2Phase()

        calculated_state = state()  # It stores new objects into the state class.
        calculated_state.t = self.t # Current temp
        calculated_state.p = self.p # Current Pressure
        calculated_state.v = self.v # Current volume
        calculated_state.u = self.u # Current internal energy
        calculated_state.h = self.h # Current enthalpy
        calculated_state.s = self.s # Current entropy
        calculated_state.x = self.x # Current quality



        self._lbl_StateProperties_1.setText(self.stProps)
        return calculated_state # returns the populated state

    def calculateProperties2(self):
        """
        Calculates the thermodynamic state variables based on specified values.
        I have thermodynamic variables:  P, T, v, h, u, s and x (7 things) from which I am choosing two.
        Possible number of permutations:  7!/5! =42.
        But, order of the two things does not matter, so 42/2=21
        PT, Pv, Ph, Pu, Ps, Px (6)
        Tv, Th, Tu, Ts, Tx (5)
        vh, vu, vs, vx (4)
        hu, hs, hx (3)
        us, ux (2)
        sx (1)
        Total of 21 cases to deal with.  I will attack them in the order shown above. This is the same as the function
        before but changes to account for the state 2.
        :return: nothing
        """
        # Step 1: read which properties are being specified from the combo boxes
        SP=[self._cmb_Property1_2.currentText()[-2:-1], self._cmb_Property2_2.currentText()[-2:-1]]
        if SP[0]==SP[1]:
            self._lbl_Warning.setText("Warning:  You cannot specify the same property twice.")
        else:
            self._lbl_Warning.setText("")
        SP[0]=SP[0].lower()
        SP[1]=SP[1].lower()

        # Step 2: select the proper case from the 21.  Note that PT is the same as TP etc.
        if SP[0]=='p' or SP[1]=='p':
            oFlipped= SP[0] != 'p'
            SP1 = SP[0] if oFlipped else SP[1]
        #case 1:  pt or tp
            if SP1=='t':
                f1=float(self._le_Property1_2.text())
                f2=float(self._le_Property2_2.text())
                self.p= f1 if not oFlipped else f2
                self.t= f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.tSat=round(self.tSat,3)  # I will compare at 3 three decimal places
                # compare T to TSat
                if self.t<self.tSat or self.t>self.tSat:
                    self.region = "sub-cooled liquid" if self.t<self.tSat else "super-heated vapor"
                    self.makeLabel_1Phase_state2()
                else:  #this is ambiguous since at saturated temperature
                    self.region = "two-phase"
                    self.stProps ="Region = {:}".format(self.region)
                    self.stProps += "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
                    self.stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat, self.t_Units)
                    self.stProps += "\nInternal Energy = uf={:0.3f}, ug={:0.3f} ({:})".format(self.uf, self.ug, self.u_Units)
                    self.stProps += "\nEnthalpy = hf={:0.3f}, hg={:0.3f} ({:})".format(self.hf,self.hg, self.h_Units)
                    self.stProps += "\nEntropy = sf={:0.3f}, sg={:0.3f} ({:})".format(self.sf, self.sg, self.s_Units)
                    self.stProps += "\nSpecific Volume = vf={:0.3f}, vg={:0.3f} ({:})".format(self.vf, self.vg, self.v_Units)
                    self.stProps += "\nQuality = unknown"
        #case 2: pv or vp
            elif SP1=='v':
                f1=float(self._le_Property1_2.text())
                f2=float(self._le_Property2_2.text())
                self.p= f1 if not oFlipped else f2
                self.v= f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.vf=round(self.vf,5)
                self.vg=round(self.vg,3)
                # compare v to vf and vg
                if self.v<self.vf or self.v>self.vg:
                    self.region = "sub-cooled liquid" if self.v<self.vf else "super-heated vapor"
                    #since I can't find properties using v, I will use fsolve to find T
                    dt=1.0 if self.v>self.vg else -1.0
                    fn1 = lambda T: self.v - self.steamTable.v_pt(self.p, T)
                    self.t = fsolve(fn1, [self.tSat + dt])[0]
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.v-self.vf)/(self.vgf)
                    self.t=self.tSat
                    self.makeLabel_2Phase_state2()
        #case 3 pu or up
            elif SP1=='u':
                f1=float(self._le_Property1_2.text())
                f2=float(self._le_Property2_2.text())
                self.p= f1 if not oFlipped else f2;
                self.u= f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u<self.uf or self.u>self.ug:
                    self.region = "sub-cooled liquid" if self.u<self.uf else "super-heated vapor"
                    #since I can't find properties using u, I will use fsolve to find T
                    dt=1.0 if self.u>self.ug else -1.0
                    fn3 = lambda T: self.u - self.steamTable.u_pt(self.p, T)
                    self.t = fsolve(fn3, [self.tSat + dt])[0]
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.u-self.uf)/(self.ugf)
                    self.t=self.tSat
                    self.makeLabel_2Phase_state2()
        #case 4 ph or hp
            elif SP1 == 'h':
                f1=float(self._le_Property1_2.text())
                f2=float(self._le_Property2_2.text())
                self.p= f1 if not oFlipped else f2;
                self.h= f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare h to hf and hg
                if self.h < self.hf or self.h > self.hg:
                    self.region = "sub-cooled liquid" if self.h < self.hf else "super-heated vapor"
                    self.t=self.steamTable.t_ph(self.p, self.h)
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.h - self.hf) / (self.hgf)
                    self.t=self.tSat
                    self.makeLabel_2Phase_state2()
        #case 5 ps or sp
            elif SP1 == 's':
                f1=float(self._le_Property1_2.text())
                f2=float(self._le_Property2_2.text())
                self.p= f1 if not oFlipped else f2;
                self.s= f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    self.t=self.steamTable.t_ps(self.p, self.s)
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase_state2()
        #case 6 px or xp
            elif SP1 == 'x':
                self.region="two-phase"
                f1=float(self._le_Property1_2.text())
                f2=float(self._le_Property2_2.text())
                self.p= f1 if not oFlipped else f2
                self.x= f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.t=self.tSat
                self.x=self.clamp(self.x,0.0,1.0)
                self.makeLabel_2Phase_state2()
        elif SP[0]=='t' or SP[1]=='t':
            oFlipped = SP[0] != 't'
            SP1 = SP[0] if oFlipped else SP[1]
        #case 7:  tv or vt
            if SP1=='v':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.v = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                self.vf=round(self.vf,5)
                self.vg=round(self.vg,3)
                # compare v to vf and vg
                if self.v<self.vf or self.v>self.vg:
                    self.region = "sub-cooled liquid" if self.v<self.vf else "super-heated vapor"
                    #since I can't find properties using v, I will use fsolve to find P
                    dp=-0.1 if self.v>self.vg else 0.1
                    fn3 = lambda P: [self.v - self.steamTable.v_pt(P, self.t)]
                    self.p = fsolve(fn3, [self.pSat + dp])[0]
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.v-self.vf)/(self.vgf)
                    self.p=self.pSat
                    self.makeLabel_2Phase_state2()
        #case 8:  tu or ut
            elif SP1=='u':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u<self.uf else "super-heated vapor"
                    #since I can't find properties using u, I will use fsolve to find P
                    dp=0.1 if self.u>self.ug else -0.1
                    fn8 = lambda P: self.u - self.steamTable.u_pt(self.t, P)
                    self.p = fsolve(fn8, [self.pSat + dp])[0]
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.u-self.uf)/(self.ugf)
                    self.p=self.pSat
                    self.makeLabel_2Phase_state2()
        #case 9:  th or ht
            elif SP1 == 'h':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.h = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare h to hf and hg
                if self.h < self.hf or self.h > self.hg:
                    self.region = "sub-cooled liquid" if self.h < self.hf else "super-heated vapor"
                    self.p=self.steamTable.p_th(self.t, self.h)
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p=self.pSat
                    self.x = (self.h - self.hf) / (self.hgf)
                    self.makeLabel_2Phase_state2()
        #case 10:  ts or st
            elif SP1 == 's':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    self.p=self.steamTable.p_ts(self.t, self.s)
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p=self.pSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase_state2()
        #case 11:  tx or xt
            elif SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region="two-phase"
                self.getSatProps_t(self.t)
                self.p=self.pSat
                self.x = float(self._le_Property2_2.text())
                self.x=self.clamp(self.x,0.0,1.0)
                self.makeLabel_2Phase_state2()
        elif SP[0] == 'v' or SP[1] == 'v':
            oFlipped = SP[0] != 'v'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 12:  vh or hv
            if SP1 == 'h':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.h = f2 if not oFlipped else f1
                def fn12(P):
                    # could be single phase or two-phase, but both v&h have to match at same x
                    self.getSatProps_p(P)
                    if self.between(self.h,self.hf, self.hg):
                        self.x=(self.h-self.hf)/self.hgf
                        return self.v-(self.vf+self.x*self.vgf)
                    # could be single phase
                    return self.v-self.steamTable.v_ph(P,self.h)
                self.p=fsolve(fn12,[1.0])[0]
                self.t=self.steamTable.t_ph(self.p, self.h)
                self.getSatProps_p(self.p)
                self.vf = round(self.vf, 5)
                self.vg = round(self.vg, 3)
                # compare v to vf and vg
                if self.v < self.vf or self.v > self.vg:
                    self.region = "sub-cooled liquid" if self.v < self.vf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.v - self.vf) / (self.vgf)
                    self.makeLabel_2Phase_state2()
            # case 13:  vu or uv
            elif SP1 == 'u':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                # use fsolve to fing P&T at this v & u
                def fn13(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.u,self.uf, self.ug):
                        self.t=self.tSat
                        self.x=(self.u-self.uf)/self.ugf
                        return [self.v-(self.vf+self.x*self.vgf),0]
                    return [self.v-self.steamTable.v_pt(PT[0],PT[1]),self.u-self.steamTable.u_pt(PT[0],PT[1])]
                props=fsolve(fn13, [1,100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.p = self.pSat
                    self.t=self.tSat
                    self.makeLabel_2Phase_state2()
            # case 14:  vs os sv
            elif SP1 == 's':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                def fn14(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x=(self.s-self.sf)/self.sgf
                        return [self.v-self.vf-self.x*self.vgf, 0.0]
                    return [self.v - self.steamTable.v_pt(PT[0], PT[1]),
                                  self.s - self.steamTable.s_pt(PT[0], PT[1])]
                props = fsolve(fn14, [1, 100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase_state2()
            # case 15:  vx or xv
            elif SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn15(p):
                    self.getSatProps_p(p)
                    return self.v -(self.vf+self.x*self.vgf)
                props=fsolve(fn15, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_2.text()),0.0,1.0)
                self.makeLabel_2Phase_state2()
        elif SP[0] == 'h' or SP[1] == 'h':
            oFlipped = SP[0] != 'h'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 16:  hu or uh
            if SP1 == 'u':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.h = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                # use fsolve to fing P&T at this v & u
                def fn16(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.u, self.uf, self.ug):
                        self.x=(self.u-self.uf)/self.ugf
                        return [self.h-self.hf-self.x*self.hgf, 0.0]
                    return [self.h-self.steamTable.h_pt(PT[0],PT[1]),self.u-self.steamTable.u_pt(PT[0],PT[1])]
                props=fsolve(fn16, [1,100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.p = self.pSat
                    self.t=self.tSat
                    self.makeLabel_2Phase_state2()
            # case 17:  hs or sh
            elif SP1 == 's':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.h = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                def fn17(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x=(self.s-self.sf)/self.sgf
                        return [self.h-self.hf-self.x*self.hgf, 0.0]
                    return [self.h-self.steamTable.h_pt(PT[0],PT[1]),self.s-self.steamTable.s_pt(PT[0],PT[1])]
                props = fsolve(fn17, [1, 100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase_state2()
            # case 18:  hx or xh
            elif SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn18(p):
                    self.getSatProps_p(p)
                    return self.h -(self.hf+self.x*self.hgf)
                props=fsolve(fn18, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_2.text()),0.0,1.0)
                self.makeLabel_2Phase_state2()
        elif SP[0] == 'u' or SP[1] == 'u':
            oFlipped = SP[0] != 'u'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 19:  us or su
            if SP1 == 's':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.u = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                def fn19(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x = (self.s - self.sf) / self.sgf
                        return [self.u - self.uf - self.x * self.ugf, 0.0]
                    return [self.u - self.steamTable.u_pt(PT[0], PT[1]),
                            self.s - self.steamTable.s_pt(PT[0], PT[1])]
                props = fsolve(fn19, [1, 100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase_state2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t=self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase_state2()
            # case 20:  ux or xu
            elif SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.u = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn20(p):
                    self.getSatProps_p(p)
                    return self.h -(self.hf+self.x*self.hgf)
                props=fsolve(fn20, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_2.text()),0.0,1.0)
                self.makeLabel_2Phase_state2()
        elif SP[0] == 's' or SP[1] == 's':
            oFlipped = SP[0] != 's'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 21:  sx or xs
            if SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.s = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn21(p):
                    self.getSatProps_p(p)
                    return self.h -(self.hf+self.x*self.hgf)
                props=fsolve(fn21, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_2.text()),0.0,1.0)
                self.makeLabel_2Phase_state2()

        calculated_state = state2()  # It stores values again in the state class
        calculated_state.t = self.t
        calculated_state.p = self.p
        calculated_state.v = self.v
        calculated_state.u = self.u
        calculated_state.h = self.h
        calculated_state.s = self.s
        calculated_state.x = self.x

        self._lbl_StateProperties_2.setText(self.stProps)
        return calculated_state #returns current values state 2

    def calculateDelta(self, state1, state2):
        """
        Calculates the difference in thermodynamic properties between two states.
        :param state1:The initial state object
        :param state2: The final state object.
        This function computes the differences in temperature, pressure, specific volume,
        internal energy, enthalpy, entropy, and quality between the two provided state objects.
        It also updates the UI to display these differences.
        :return: state: A state object containing the differences in properties between state1 and state2
        """
        # Calculate the raw differences
        delta_state = state()
        delta_state.t = state2.t - state1.t
        delta_state.p = state2.p - state1.p
        delta_state.v = state2.v - state1.v
        delta_state.u = state2.u - state1.u
        delta_state.h = state2.h - state1.h
        delta_state.s = state2.s - state1.s
        delta_state.x = state2.x - state1.x


        # Prepare the text with units and deltas
        deltaText = ""
        deltaText += f"ΔPressure = {delta_state.p:.3f} {self.p_Units}\n"
        deltaText += f"ΔTemperature = {delta_state.t:.3f} {self.t_Units}\n"
        deltaText += f"ΔInternal Energy = {delta_state.u:.3f} {self.u_Units}\n"
        deltaText += f"ΔEnthalpy = {delta_state.h:.3f} {self.h_Units}\n"
        deltaText += f"ΔEntropy = {delta_state.s:.3f} {self.s_Units}\n"
        deltaText += f"ΔSpecific Volume = {delta_state.v:.5f} {self.v_Units}\n"
        deltaText += f"ΔQuality = {delta_state.x:.3f}" if delta_state.x is not None else "ΔQuality = unknown"

        # Set the label text
        self._lbl_StateProperties_State_change.setText(deltaText)

        return delta_state


#endregion

#region function definitions
def main():
    """
    There seems to be a bug, that makes it so that to properly run the program, the first time
    you must click calculate with the standard values. Then you can change it works. But sometimes some values cause
    program to close, I believe it may be from bad values but there might be a small error in the code that I can't find.

    Entry point of the application that initializes and starts the PyQt application loop.

    This function performs the following actions to launch the thermodynamic state calculator application:

    - Checks if a QApplication instance already exists to avoid creating a second instance.
    - Creates a new QApplication instance if one does not already exist, passing in any command-line arguments.
    - Connects the application's 'aboutToQuit' signal to the 'deleteLater' slot to ensure a clean exit.
    - Creates an instance of the main window class, which sets up the UI and its functionality.
    - Starts the event loop of the application, which waits for user interactions such as button clicks,
      and processes these events as defined in the application's slots and signals.
    - Exits the application once the event loop is terminated and returns the exit status to the system.

    This main function essentially bootstraps the application, setting up the necessary components
    for the graphical user interface and entering the main event loop that handles user interactions.
    :return:
    """
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    sys.exit(app.exec_())
    pass
#end region

#region function calls
if __name__=="__main__":
    main()