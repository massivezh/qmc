#!/usr/bin/env python

# Simple GUI for qmc.py
# FIXME Experimental - doesn't do any check on what is passed as input ;)
#
#  Copyright (C) 2011 Marcello Pogliani
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import qmc
import sys
import signal
from PyQt4.Qt import *

signal.signal(signal.SIGINT, signal.SIG_DFL)

# TODO refactor the library to allow GUI to output some intermediate steps!

if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = QWidget()
    
    widget.resize(450, 350)
    widget.setWindowTitle('Quine McCluskey Algorithm')
    
    layout = QGridLayout(widget)
    widget.setLayout(layout)
    
    # widgets
    go = QPushButton('GO!', widget)
    reset = QPushButton('Reset', widget)
    add_function = QPushButton('Add function', widget)
    costfun_selector = QButtonGroup(widget)
    costfun_selector_literals = QRadioButton('# of literals', widget)
    costfun_selector_implicants = QRadioButton('# of implicants', widget)
    costfun_selector.addButton(costfun_selector_literals)
    costfun_selector.addButton(costfun_selector_implicants)
    costfun_selector_literals.setChecked(True) # default cost function
    cost = QLCDNumber(widget)
    result = QTextEdit(widget)
    insert_pane = QTableWidget(1, 2, widget);
    insert_pane.setHorizontalHeaderLabels(['ONset', 'DCset'])
    
    # bind widgets to layout
    layout.addWidget (insert_pane, 1, 1, 1, 4)
    
    layout.addWidget(add_function, 2, 1, 1, 1)
    layout.addWidget(go, 2, 2, 1, 2)
    layout.addWidget(reset, 2, 4, 1, 1)
    
    layout.addWidget(QLabel('Cost function:', widget), 3, 1, 1, 2)
    layout.addWidget(costfun_selector_implicants, 4, 1, 1, 2)
    layout.addWidget(costfun_selector_literals, 5, 1, 1, 2)
    layout.addWidget(QLabel('Computed cost:', widget), 6, 1, 2, 1)
    layout.addWidget(cost, 6, 2, 2, 1)
    
    layout.addWidget(result, 3, 3, 5, 2)
    
    def addFunction():
        insert_pane.setRowCount(insert_pane.rowCount()+1)
        
    def toList(obj):
        if obj == None:
            l = []
        else:
            s = obj.text().toAscii()
            l = s.split(',')
            l = [i.toInt()[0] for i in l]
        return l
    
    def startMinimization():
        lof = []
        for row in range(insert_pane.rowCount()):
            curf_onset = toList(insert_pane.item(row, 0))
            curf_dcset = toList(insert_pane.item(row, 1))
            if curf_onset != []:
                lof.append(qmc.QmcFunction(curf_onset, curf_dcset))
        if costfun_selector_literals.isChecked():
            costf = qmc.LITERALS_COST_FUNCTION
        elif costfun_selector_implicants.isChecked():
            costf = qmc.IMPLICANTS_COST_FUNCTION
        if lof != []:
            qmc.VERBOSE = False # no debug printfs when running from the GUI!
            q = qmc.QuineMcCluskey(lof, costf)
            q.findPrimeImplicants()
            q.simplify()
            result.setText(str(q.sol))
            cost.display(q.sol.getCost())
        else:
            result.setText("Input is empty!")
    
    def clearAll():
        insert_pane.setRowCount(1)
        insert_pane.clearContents()
        result.clear()
        cost.display(0)
        pass
    
    widget.connect(add_function, SIGNAL('clicked()'), addFunction)
    widget.connect(go, SIGNAL('clicked()'), startMinimization)
    widget.connect(reset, SIGNAL('clicked()'), clearAll)
    
    widget.show()
    sys.exit(app.exec_())
