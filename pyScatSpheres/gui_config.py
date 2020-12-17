import easygui

url = 'http://localhost:8002/scattering/linear_array_qdot_sphere/'
help_msg='''
                ###### hard sphere scattering ########
This is a gui for the scattering of a linear array of spheres with constant potential.
More info on the formulation of the problem at :
%s

It is running with matplotlib, so be gentle with it.
Use the SPACEBAR to refresh the view if something goes wrong especially after resizing the window.

## description
In this default setup the problem has been solved for N=2 spheres at few normalized radius ka
and normalized distance kd.

##Graph description :
top-left     : The total scattering cross section
top-right    : The scattering amplitude/differential cross section.The checkboxes on the left side allow you to select the real,imaginary part,phase,magnitude or magnitude square.
bottom-left  : The spherical expansion coefficients in the complex plane for the selected value of
bottom-right : the spheres setup. click on the blue button to display the near field

## GUI usage:''' %url

help_cmd = '''
Click on a (ka,kd) point on the scattering cross section to show the data corresponding to that point
    keyboard commands :
arrows    : move the active selected (ka,kd) horizontally or vertically
enter     : change settings
F1        : show this help
F5        : refresh view
space     : refresh
'''

def multenterbox(msg,title,fieldValues,fieldNames):
    fieldValues = easygui.multenterbox(msg, title, fieldNames,fieldValues)
    while True:
        if fieldValues is None:
            break
        errs = list()
        for n, v in zip(fieldNames, fieldValues):
            if v.strip() == "":errs.append('"{}" is a required field.'.format(n))
        if not len(errs):
            break
        fieldValues = easygui.multenterbox("\n".join(errs), title, fieldNames, fieldValues)
    return fieldValues
