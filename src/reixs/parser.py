import re
import parser
from .util import doesMatchPattern
from numpy import log as ln
from numpy import log10 as log
from numpy import exp


def math_stream(formula, data, arg, get_data, XAS_streams=None, is_XAS=False, background=None, REIXSObj=None):
    # Split the user input string at all mathematical operations
    # Allow "( ) * / + -" as math

    pattern = '[\(+\-*^/\)]'
    split_expr = re.split(pattern, formula)

    # Place string literals in dict if cannot be converted to float
    # drop all empty strings from re.split
    quantity_str_dict = dict()

    for i, string in enumerate(split_expr):
        if string != "":
            try:
                float(string)
            except:
                # Use math expressions to allow logs and exps
                math_expressions = ['ln', 'log', 'exp']

                if string in math_expressions:
                    pass

                else:
                    # Assign generic "val{i}" key to string literal in compliance with
                    # python supported syntax for variables
                    quantity_str_dict[f"val{i}"] = string

    # Parser does not support special string literals (ROIs) due to inproper python variable naming syntax
    # Replace them with generic key as per dictionary
    # Create local variable and assign corresponding data -- needed as eval interprets all input as variable
    # This is indeed local to this function only and cannot be accessed from loadSCAscans
    for k, v in quantity_str_dict.items():
        formula = formula.replace(v, k)

        # Ensure that in XASLoader all quantities are normalized by mesh
        # per definition of what XAS is
        # Exclude those quantities from normalization that have been normalized elsewhere

        if is_XAS == False:
            locals()[k] = get_data(v, data, arg, background, REIXSObj)
        else:
            if doesMatchPattern(v, XAS_streams):
                locals()[k] = get_data(v, data, arg, background, REIXSObj)
            else:
                numerator = get_data(v, data, arg, background, REIXSObj)
                mesh = data[arg].mesh_current
                locals()[k] = numerator/mesh

    # Return the calculated result
    code = parser.expr(formula).compile()
    return eval(code)
