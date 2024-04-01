# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ThermoStateCalc_update.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui__frm_StateCalculator(object):
    def setupUi(self, _frm_StateCalculator):
        """
        This is the transformed version of the Ui created using Qt designer, It corresponds to the problem from the
        homework
        :param _frm_StateCalculator:
        :return: None
        """
        _frm_StateCalculator.setObjectName("_frm_StateCalculator")
        _frm_StateCalculator.resize(747, 418)
        self.verticalLayout = QtWidgets.QVBoxLayout(_frm_StateCalculator)
        self.verticalLayout.setObjectName("verticalLayout")
        self._grp_Units = QtWidgets.QGroupBox(_frm_StateCalculator)
        self._grp_Units.setObjectName("_grp_Units")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self._grp_Units)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self._rdo_SI = QtWidgets.QRadioButton(self._grp_Units)
        self._rdo_SI.setChecked(True)
        self._rdo_SI.setObjectName("_rdo_SI")
        self.horizontalLayout_2.addWidget(self._rdo_SI)
        self._rdo_English = QtWidgets.QRadioButton(self._grp_Units)
        self._rdo_English.setObjectName("_rdo_English")
        self.horizontalLayout_2.addWidget(self._rdo_English)
        self.verticalLayout.addWidget(self._grp_Units)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self._grp_SpecifiedProperties = QtWidgets.QGroupBox(_frm_StateCalculator)
        self._grp_SpecifiedProperties.setObjectName("_grp_SpecifiedProperties")
        self.gridLayout = QtWidgets.QGridLayout(self._grp_SpecifiedProperties)
        self.gridLayout.setObjectName("gridLayout")
        self.State1 = QtWidgets.QGroupBox(self._grp_SpecifiedProperties)
        self.State1.setObjectName("State1")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.State1)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self._lbl_Property1_1 = QtWidgets.QLabel(self.State1)
        self._lbl_Property1_1.setObjectName("_lbl_Property1_1")
        self.gridLayout_2.addWidget(self._lbl_Property1_1, 0, 0, 1, 1)
        self._lbl_Property2_1 = QtWidgets.QLabel(self.State1)
        self._lbl_Property2_1.setObjectName("_lbl_Property2_1")
        self.gridLayout_2.addWidget(self._lbl_Property2_1, 0, 2, 1, 1)
        self._cmb_Property1_1 = QtWidgets.QComboBox(self.State1)
        self._cmb_Property1_1.setObjectName("_cmb_Property1_1")
        self._cmb_Property1_1.addItem("")
        self._cmb_Property1_1.addItem("")
        self._cmb_Property1_1.addItem("")
        self._cmb_Property1_1.addItem("")
        self._cmb_Property1_1.addItem("")
        self._cmb_Property1_1.addItem("")
        self._cmb_Property1_1.addItem("")
        self.gridLayout_2.addWidget(self._cmb_Property1_1, 1, 0, 1, 2)
        self._cmb_Property2_1 = QtWidgets.QComboBox(self.State1)
        self._cmb_Property2_1.setObjectName("_cmb_Property2_1")
        self._cmb_Property2_1.addItem("")
        self._cmb_Property2_1.addItem("")
        self._cmb_Property2_1.addItem("")
        self._cmb_Property2_1.addItem("")
        self._cmb_Property2_1.addItem("")
        self._cmb_Property2_1.addItem("")
        self._cmb_Property2_1.addItem("")
        self.gridLayout_2.addWidget(self._cmb_Property2_1, 1, 2, 1, 2)
        self._le_Property1_1 = QtWidgets.QLineEdit(self.State1)
        self._le_Property1_1.setObjectName("_le_Property1_1")
        self.gridLayout_2.addWidget(self._le_Property1_1, 2, 0, 1, 1)
        self._lbl_Property1_1_Units = QtWidgets.QLabel(self.State1)
        self._lbl_Property1_1_Units.setObjectName("_lbl_Property1_1_Units")
        self.gridLayout_2.addWidget(self._lbl_Property1_1_Units, 2, 1, 1, 1)
        self._le_Property2_1 = QtWidgets.QLineEdit(self.State1)
        self._le_Property2_1.setObjectName("_le_Property2_1")
        self.gridLayout_2.addWidget(self._le_Property2_1, 2, 2, 1, 1)
        self._lbl_Property2_1_Units = QtWidgets.QLabel(self.State1)
        self._lbl_Property2_1_Units.setObjectName("_lbl_Property2_1_Units")
        self.gridLayout_2.addWidget(self._lbl_Property2_1_Units, 2, 3, 1, 1)
        self.gridLayout.addWidget(self.State1, 0, 0, 2, 1)
        self._pb_Calculate = QtWidgets.QPushButton(self._grp_SpecifiedProperties)
        self._pb_Calculate.setObjectName("_pb_Calculate")
        self.gridLayout.addWidget(self._pb_Calculate, 3, 0, 1, 2)
        self._lbl_Warning = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Warning.setText("")
        self._lbl_Warning.setObjectName("_lbl_Warning")
        self.gridLayout.addWidget(self._lbl_Warning, 2, 0, 1, 2)
        self.State2 = QtWidgets.QGroupBox(self._grp_SpecifiedProperties)
        self.State2.setObjectName("State2")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.State2)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self._lbl_Property2_2_Units = QtWidgets.QLabel(self.State2)
        self._lbl_Property2_2_Units.setObjectName("_lbl_Property2_2_Units")
        self.gridLayout_3.addWidget(self._lbl_Property2_2_Units, 2, 3, 1, 1)
        self._lbl_Property1_2 = QtWidgets.QLabel(self.State2)
        self._lbl_Property1_2.setObjectName("_lbl_Property1_2")
        self.gridLayout_3.addWidget(self._lbl_Property1_2, 0, 0, 1, 1)
        self._lbl_Property1_2_Units = QtWidgets.QLabel(self.State2)
        self._lbl_Property1_2_Units.setObjectName("_lbl_Property1_2_Units")
        self.gridLayout_3.addWidget(self._lbl_Property1_2_Units, 2, 1, 1, 1)
        self._lbl_Property2_2 = QtWidgets.QLabel(self.State2)
        self._lbl_Property2_2.setObjectName("_lbl_Property2_2")
        self.gridLayout_3.addWidget(self._lbl_Property2_2, 0, 2, 1, 1)
        self._cmb_Property1_2 = QtWidgets.QComboBox(self.State2)
        self._cmb_Property1_2.setObjectName("_cmb_Property1_2")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self.gridLayout_3.addWidget(self._cmb_Property1_2, 1, 0, 1, 2)
        self._le_Property2_2 = QtWidgets.QLineEdit(self.State2)
        self._le_Property2_2.setObjectName("_le_Property2_2")
        self.gridLayout_3.addWidget(self._le_Property2_2, 2, 2, 1, 1)
        self._cmb_Property2_2 = QtWidgets.QComboBox(self.State2)
        self._cmb_Property2_2.setObjectName("_cmb_Property2_2")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self.gridLayout_3.addWidget(self._cmb_Property2_2, 1, 2, 1, 2)
        self._le_Property1_2 = QtWidgets.QLineEdit(self.State2)
        self._le_Property1_2.setObjectName("_le_Property1_2")
        self.gridLayout_3.addWidget(self._le_Property1_2, 2, 0, 1, 1)
        self.gridLayout.addWidget(self.State2, 0, 1, 2, 1)
        self.verticalLayout.addWidget(self._grp_SpecifiedProperties)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self._grp_StateProperties = QtWidgets.QGroupBox(_frm_StateCalculator)
        self._grp_StateProperties.setObjectName("_grp_StateProperties")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self._grp_StateProperties)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self._lbl_SatVapProps = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_SatVapProps.setText("")
        self._lbl_SatVapProps.setObjectName("_lbl_SatVapProps")
        self.horizontalLayout.addWidget(self._lbl_SatVapProps)
        self.State_1 = QtWidgets.QGroupBox(self._grp_StateProperties)
        self.State_1.setObjectName("State_1")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.State_1)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self._lbl_State_1 = QtWidgets.QLabel(self.State_1)
        self._lbl_State_1.setObjectName("_lbl_State_1")
        self.verticalLayout_2.addWidget(self._lbl_State_1)
        self._lbl_StateProperties_1 = QtWidgets.QLabel(self.State_1)
        self._lbl_StateProperties_1.setObjectName("_lbl_StateProperties_1")
        self.verticalLayout_2.addWidget(self._lbl_StateProperties_1)
        self.horizontalLayout.addWidget(self.State_1)
        self.State_2 = QtWidgets.QGroupBox(self._grp_StateProperties)
        self.State_2.setObjectName("State_2")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.State_2)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self._lbl_State_2 = QtWidgets.QLabel(self.State_2)
        self._lbl_State_2.setObjectName("_lbl_State_2")
        self.verticalLayout_3.addWidget(self._lbl_State_2)
        self._lbl_StateProperties_2 = QtWidgets.QLabel(self.State_2)
        self._lbl_StateProperties_2.setObjectName("_lbl_StateProperties_2")
        self.verticalLayout_3.addWidget(self._lbl_StateProperties_2)
        self.horizontalLayout.addWidget(self.State_2)
        self.State_Change = QtWidgets.QGroupBox(self._grp_StateProperties)
        self.State_Change.setObjectName("State_Change")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.State_Change)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self._lbl_StateProperties_State_change = QtWidgets.QLabel(self.State_Change)
        self._lbl_StateProperties_State_change.setObjectName("_lbl_StateProperties_State_change")
        self.verticalLayout_4.addWidget(self._lbl_StateProperties_State_change)
        self.horizontalLayout.addWidget(self.State_Change)
        self._lbl_SatLiqProps = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_SatLiqProps.setText("")
        self._lbl_SatLiqProps.setObjectName("_lbl_SatLiqProps")
        self.horizontalLayout.addWidget(self._lbl_SatLiqProps)
        self.verticalLayout.addWidget(self._grp_StateProperties)

        self.retranslateUi(_frm_StateCalculator)
        self._cmb_Property2_1.setCurrentIndex(1)
        self._cmb_Property2_2.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(_frm_StateCalculator)

    def retranslateUi(self, _frm_StateCalculator):
        _translate = QtCore.QCoreApplication.translate
        _frm_StateCalculator.setWindowTitle(_translate("_frm_StateCalculator", "Thermodynamic State Calculator"))
        self._grp_Units.setTitle(_translate("_frm_StateCalculator", "System of Units"))
        self._rdo_SI.setText(_translate("_frm_StateCalculator", "SI"))
        self._rdo_English.setText(_translate("_frm_StateCalculator", "English"))
        self._grp_SpecifiedProperties.setTitle(_translate("_frm_StateCalculator", "Specified Properties"))
        self.State1.setTitle(_translate("_frm_StateCalculator", "State 1"))
        self._lbl_Property1_1.setText(_translate("_frm_StateCalculator", "Property 1"))
        self._lbl_Property2_1.setText(_translate("_frm_StateCalculator", "Property 2"))
        self._cmb_Property1_1.setItemText(0, _translate("_frm_StateCalculator", "Pressure (p)"))
        self._cmb_Property1_1.setItemText(1, _translate("_frm_StateCalculator", "Temperature (T)"))
        self._cmb_Property1_1.setItemText(2, _translate("_frm_StateCalculator", "Quality (x)"))
        self._cmb_Property1_1.setItemText(3, _translate("_frm_StateCalculator", "Specific Internal Energy (u)"))
        self._cmb_Property1_1.setItemText(4, _translate("_frm_StateCalculator", "Specific Enthalpy (h)"))
        self._cmb_Property1_1.setItemText(5, _translate("_frm_StateCalculator", "Specific Volume (v)"))
        self._cmb_Property1_1.setItemText(6, _translate("_frm_StateCalculator", "Specific Entropy (s)"))
        self._cmb_Property2_1.setItemText(0, _translate("_frm_StateCalculator", "Pressure (p)"))
        self._cmb_Property2_1.setItemText(1, _translate("_frm_StateCalculator", "Temperature (T)"))
        self._cmb_Property2_1.setItemText(2, _translate("_frm_StateCalculator", "Quality (x)"))
        self._cmb_Property2_1.setItemText(3, _translate("_frm_StateCalculator", "Specific Internal Energy (u)"))
        self._cmb_Property2_1.setItemText(4, _translate("_frm_StateCalculator", "Specific Enthalpy (h)"))
        self._cmb_Property2_1.setItemText(5, _translate("_frm_StateCalculator", "Specific Volume (v)"))
        self._cmb_Property2_1.setItemText(6, _translate("_frm_StateCalculator", "Specific Entropy (s)"))
        self._le_Property1_1.setText(_translate("_frm_StateCalculator", "1.0"))
        self._lbl_Property1_1_Units.setText(_translate("_frm_StateCalculator", "Bar"))
        self._le_Property2_1.setText(_translate("_frm_StateCalculator", "100.0"))
        self._lbl_Property2_1_Units.setText(_translate("_frm_StateCalculator", "C"))
        self._pb_Calculate.setText(_translate("_frm_StateCalculator", "Calculate"))
        self.State2.setTitle(_translate("_frm_StateCalculator", "State 2"))
        self._lbl_Property2_2_Units.setText(_translate("_frm_StateCalculator", "C"))
        self._lbl_Property1_2.setText(_translate("_frm_StateCalculator", "Property 1"))
        self._lbl_Property1_2_Units.setText(_translate("_frm_StateCalculator", "Bar"))
        self._lbl_Property2_2.setText(_translate("_frm_StateCalculator", "Property 2"))
        self._cmb_Property1_2.setItemText(0, _translate("_frm_StateCalculator", "Pressure (p)"))
        self._cmb_Property1_2.setItemText(1, _translate("_frm_StateCalculator", "Temperature (T)"))
        self._cmb_Property1_2.setItemText(2, _translate("_frm_StateCalculator", "Quality (x)"))
        self._cmb_Property1_2.setItemText(3, _translate("_frm_StateCalculator", "Specific Internal Energy (u)"))
        self._cmb_Property1_2.setItemText(4, _translate("_frm_StateCalculator", "Specific Enthalpy (h)"))
        self._cmb_Property1_2.setItemText(5, _translate("_frm_StateCalculator", "Specific Volume (v)"))
        self._cmb_Property1_2.setItemText(6, _translate("_frm_StateCalculator", "Specific Entropy (s)"))
        self._le_Property2_2.setText(_translate("_frm_StateCalculator", "100.0"))
        self._cmb_Property2_2.setItemText(0, _translate("_frm_StateCalculator", "Pressure (p)"))
        self._cmb_Property2_2.setItemText(1, _translate("_frm_StateCalculator", "Temperature (T)"))
        self._cmb_Property2_2.setItemText(2, _translate("_frm_StateCalculator", "Quality (x)"))
        self._cmb_Property2_2.setItemText(3, _translate("_frm_StateCalculator", "Specific Internal Energy (u)"))
        self._cmb_Property2_2.setItemText(4, _translate("_frm_StateCalculator", "Specific Enthalpy (h)"))
        self._cmb_Property2_2.setItemText(5, _translate("_frm_StateCalculator", "Specific Volume (v)"))
        self._cmb_Property2_2.setItemText(6, _translate("_frm_StateCalculator", "Specific Entropy (s)"))
        self._le_Property1_2.setText(_translate("_frm_StateCalculator", "1.0"))
        self._grp_StateProperties.setTitle(_translate("_frm_StateCalculator", "State Properties"))
        self.State_1.setTitle(_translate("_frm_StateCalculator", "State 1"))
        self._lbl_State_1.setText(_translate("_frm_StateCalculator", "State:  saturated"))
        self._lbl_StateProperties_1.setText(_translate("_frm_StateCalculator", "Pressure = 1000 kPa\n"
"Temperature = 100 C\n"
"X = 1.0"))
        self.State_2.setTitle(_translate("_frm_StateCalculator", "State 2"))
        self._lbl_State_2.setText(_translate("_frm_StateCalculator", "State:  saturated"))
        self._lbl_StateProperties_2.setText(_translate("_frm_StateCalculator", "Pressure = 1000 kPa\n"
"Temperature = 100 C\n"
"X = 1.0"))
        self.State_Change.setTitle(_translate("_frm_StateCalculator", "State Change"))
        self._lbl_StateProperties_State_change.setText(_translate("_frm_StateCalculator", "Pressurechange = 1000 kPa\n"
"Temperaturechange = 100 C\n"
"Xchange = 1.0"))