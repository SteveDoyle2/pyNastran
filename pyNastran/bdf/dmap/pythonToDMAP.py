import re


def get_comment(line):
    sline = line.rstrip().split('#')
    cmd = sline[0].strip()
    if len(sline) == 1:
        comment = ''
    elif len(sline) == 2:
        comment = '$ ' + sline[1]
    else:
        raise NotImplementedError(line)

    nspaces = len(sline[0]) - len(cmd)
    return (nspaces, cmd, comment)

from builtins import INT, FLOAT, DOUBLE, STRING, ZEROS, ONES, ARRAY


class PythonToDMAP(object):
    def __init__(self, pyName):
        self.pyName = pyName
        self.commands = []
        self.globals = {}  # start with a letter
        self.tempvars = {}  # start with a number
        self.varCounter = 0

    def run(self):
        f = open(self.pyName, 'r')
        for line in f:
            if line.strip():
                self.parse_line(line)
        f.close()

    def parse_line(self, line):
        (nspaces, cmd, comment) = get_comment(line)

        if '+=' in cmd:
            operator = '+='
            (variable, expr) = cmd.split(operator)
        elif '-=' in cmd:
            operator = '-='
            (variable, expr) = cmd.split(operator)
        elif '*=' in cmd:
            operator = '*='
            (variable, expr) = cmd.split(operator)
        elif '/=' in cmd:
            operator = '/='
            (variable, expr) = cmd.split(operator)
        elif '=' in cmd:
            operator = '='
            (variable, expr) = cmd.split(operator)
        else:
            print "cmd = |%s| not parsed" % (cmd)
            operator = None
        print "comment = ", comment
        if operator is not None:
            variable = variable.strip()
            expr = expr.strip()
            (operation, Type) = self.parse_operation(variable, operator, expr)

            self.write_operation(operation, Type, variable,
                                 operator, expr, comment)
        elif comment:
            self.commands.append(comment)
        #elif 'while' in cmd:
        #    operator =

    def write_sub_operation(self, operation, var, operator, expr):
        msg = "var=|%s| operator=|%s| expr=|%s|" % (var, operator, expr)
        if operator in ['+=', '-=', '*=', '/=']:
            operator0 = operator[0]
            call = '%s = %s %s %s' % (var, var, operator0, expr.upper())
        elif operator == '=':
            call = '%s = %s' % (var.upper(), expr.upper())
        else:
            print type(operation)
            raise NotImplementedError(msg)
        return call

    def write_operation(self, operation, Type, var, operator, expr, comment):
        if comment:
            self.commands.append('$ ' + comment)
        msg = "var=|%s| operator=|%s| expr=|%s|" % (var, operator, expr)
        if operation is None:
            msg = "var=|%s| operator=|%s| expr=|%s|" % (var, operator, expr)
            raise NotImplementedError(msg)

        assert not isinstance(Type, tuple), Type
        if isinstance(operation, INT) or isinstance(Type, INT):
            call = 'TYPE PARM,,I,N, %s $' % (var.upper())
            print call
            self.commands.append(call)
            call = self.write_sub_operation(operation, var, operator, expr)
        elif isinstance(operation, FLOAT) or isinstance(Type, FLOAT):
            call = 'TYPE PARM,,RS,N, %s $' % (var.upper())
            print call
            self.commands.append(call)
            call = self.write_sub_operation(operation, var, operator, expr)
        elif isinstance(operation, str):  # sum of integer + float
            (operation2, Type2) = self.parse_operation(var, '=', expr)
            if isinstance(Type2, STRING):
                expr2 = eval(eval(expr, self.globals))
                print "expr  = |%s|" % (expr)
                print "expr2 = |%s|" % (expr2)
                call = 'TYPE PARM,,CHAR%s,N, %s $' % (len(expr2), var.upper())
                print call
                self.commands.append(call)
            elif isinstance(Type2, DOUBLE):
                #expr2 = eval(eval(expr,self.globals))
                print "expr  = |%s|" % (expr)
                call = expr.upper()
                #print "expr2 = |%s|" %(expr2)
                call = 'TYPE PARM,,RD,N, %s $' % (var.upper())
                print call
                self.commands.append(call)
            else:
                print type(operation)
                raise NotImplementedError(msg)

            #print operation2,type(Type2)
            #asd
            #call = self.write_sub_operation(operation, var, operator, expr2)
            call = expr.upper()
        elif isinstance(operation, STRING):
            call = 'TYPE PARM,,CHAR%s,N, %s $' % (
                operation.length(self.globals) - 2, var.upper())
            print call
            self.commands.append(call)
            call = self.write_sub_operation(operation, var, operator, expr)
        elif isinstance(operation, ZEROS):
            call = 'TYPE PARM,,CHAR%s,N, %s $' % (
                operation.length(self.globals), var.upper())
            call = self.write_sub_operation(operation, var, operator, expr)
        elif isinstance(operation, DOUBLE):
            call = 'TYPE PARM,,RD,N, %s $' % (var.upper())
            print call
            self.commands.append(call)
            call = str(operation.real)
            #call = self.write_sub_operation(operation, var, operator, expr)
            print "call2 = ", call
        else:
            print type(operation)
            raise NotImplementedError(msg)
        self.commands.append(call)
        print call + '\n'

    def parse_assign(self, variable, expr):
        if expr[0] in ['+', '-']:
            sign = expr[0]
            expr = expr[1:]
        else:
            sign = ''
        sign = sign.strip('+')  # dont write extra + signs

        while '(' in expr[0]:  # get rid of extra parentheses
            expr = expr[1:-1]

        expr = expr.strip()
        if '(' in expr:  # order of operations or a function
            ### @todo not done

            i = expr.index('(')
            (pre, aft) = (expr[:i], expr[i:])
            if pre.isalpha():  # @todo function, not done
                pass
            else:  # order of operations
                self.parse_assign(self.varCounter, pre)
                self.varCounter += 1
            ###
            # stuff goes here...
            if 'array(' in expr or 'zeros(' in expr or 'ones(' in expr:
                expr_val = eval(expr.upper())
                expr_val.name = variable
                Type = expr_val
            elif 'double(' in expr:
                expr_val = eval(expr.upper())
                expr_val.name = variable
                Type = expr_val
            else:
                raise NotImplementedError('function...%s' % (expr))
            ###
            #print "expr_val = ", expr_val
        else:
            if '+' in expr:
                operator = '+'
            elif '-' in expr:
                operator = '-'
            elif '**' in expr:
                operator = '**'
            elif '*' in expr:
                operator = '*'
            elif '//' in expr:
                operator = '//'
            elif '/' in expr:
                operator = '/'
            else:
                operator = None

            if operator is not None:
                i = expr.index(operator)

                (pre, aft) = (expr[:i].strip(), expr[i:].strip())
                (expr1, Type1) = self.parse_assign(self.varCounter, pre)
                (expr2, Type2) = self.parse_assign(self.varCounter + 1, aft)
                #print "Type1=%s Type2=%s" %(Type1, Type2)
                if isinstance(Type1, INT) and isinstance(Type2, INT):
                    Type = expr1
                elif isinstance(Type1, DOUBLE) or isinstance(Type2, DOUBLE):
                    Type = expr2
                elif isinstance(Type1, FLOAT) or isinstance(Type2, FLOAT):
                    Type = expr2
                elif isinstance(Type1, STRING) and isinstance(Type2, STRING):
                    Type = expr1
                elif isinstance(Type1, ZEROS) and isinstance(Type2, ZEROS):
                    Type = expr1
                else:
                    print "Type1=%s Type2=%s" % (Type1.__class__.__name__,
                                                 Type2.__class__.__name__)
                    asf
                expr_val = expr1 + expr2
            else:
                (expr_val, Type) = self.parse_value(variable, expr)
        return (expr_val, Type)

    def parse_value(self, variable, expr):
        #print "parse_value expr = |%s|" %(expr)
        if expr.isalpha():
            expr_value = self.globals[expr]
        else:
            eval_expr = eval(expr)
            if isinstance(eval_expr, int):
                self.globals[variable] = INT(variable, expr)
            elif isinstance(eval_expr, float):
                self.globals[variable] = FLOAT(variable, expr)
            elif isinstance(eval_expr, complex):
                self.globals[variable] = COMPLEX(variable, expr)
            elif isinstance(eval_expr, str):
                self.globals[variable] = STRING(variable, expr)
            else:
                print "***eval_expr = ", eval_expr
                asdf
            ###
            expr_value = self.globals[variable]
        ###
        return (expr_value, expr_value)

    def parse_operation(self, variable, operator, expr):
        #print "var=|%s| operator=|%s| expr=|%s|" %(variable, operator, expr)
        if operator == '=':
            (expr_val, Type) = self.parse_assign(variable, expr)
            self.globals[variable] = expr_val
        return (expr_val, Type)

    def write_dmap(self, bdf):
        f = open(bdf, 'wb')
        for line in self.commands:
            f.write(line + '\n')

if __name__ == '__main__':
    comp = PythonToDMAP('spike.py')
    comp.run()
    comp.write_dmap('spike.bdf')
