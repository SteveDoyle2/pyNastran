function field = print_float_8(value)
    % Prints a float in nastran 8-character width syntax using the
    % highest precision possbile.

    if isnan(value)
        field = '        ';
        return  % field
    elseif value == 0.0
        field = sprintf('%8s', '0.');
        return  % field
    elseif value > 0.  % positive, not perfect...
        if value < 5e-8
            field = print_scientific_8(value);
            return  % field
        elseif value < 0.001
            field = print_scientific_8(value);
            field2 = sprintf("%8.7f", value);  % small value

            %field2 = field2.strip('0 ');
            field2 = regexprep(field2, '^[ 0]+|[ 0]+$', '');

            field1 = field.replace('-', 'e-');

            if field2 == '.'
                field = print_scientific_8(value);
                return  % field
            end
            if length(field2) <= 8 && float(field1) == float(field2)
                field = field2;
                %field = field.strip(' 0');
                field = regexprep(field, '^[ 0]+|[ 0]+$', '');
            end
        %elif value < 0.1
            %field = sprintf("%8.7f", value);
        elseif value < 1.
            field = sprintf("%8.7f", value);  % same as before...
        elseif value < 10.
            field = sprintf("%8.6f", value);
        elseif value < 100.
            field = sprintf("%8.5f", value);
        elseif value < 1000.
            field = sprintf("%8.4f", value);
        elseif value < 10000.
            field = sprintf("%8.3f", value);
        elseif value < 100000.
            field = sprintf("%8.2f", value);
        elseif value < 1000000.
            field = sprintf("%8.1f", value);
        else  % big value
            field = sprintf("%8.1f", value);
            %if field.index('.') < 8
            idx = strfind(field, '.');
            if idx(1) < 8
                field = sprintf('%8.1f', round(value));
                field = field(1:8);
                %assert '.' != field[0], field
            else
                field = print_scientific_8(value);
            end
            return  % field
        end
    else
        if value > -5e-7
            field = print_scientific_8(value);
            return  % field
        elseif value > -0.01  % -0.001
            field = print_scientific_8(value);
            field2 = sprintf("%8.6f", value);  % small value
            %field2 = field2.strip('0 ');
            field2 = regexprep(field2, '^[ 0]+|[ 0]+$', '');

            % get rid of the first minus sign, add it on afterwards
            % field1 = '-' + field.strip(' 0-').replace('-', 'e-');
            field1 = replace(regexprep(field, '^[ 0-]+|[ 0-]+$', ''), '-', 'e-');

            if length(field2) <= 8 && str2double(field1) == str2double(field2)
                field = field2.rstrip(' 0');
                field = field.replace('-0.', '-.');
            end

        %elseif value > -0.1
            % -0.01 >x>-0.1...should be 5 (maybe scientific...)
            %field = "%8.6f" % value
            %field = field.replace('-0.', '-.')
        elseif value > -1.
            % -0.1  >x>-1.....should be 6, but the baseline 0 is kept...
            field = sprintf("%8.6f", value);
            field = field.replace('-0.', '-.');
        elseif value > -10.
            field = sprintf("%8.5f", value);   % -1    >x>-10
        elseif value > -100.
            field = sprintf("%8.4f", value);   % -10   >x>-100
        elseif value > -1000.
            field = sprintf("%8.3f", value);   % -100  >x>-1000
        elseif value > -10000.
            field = sprintf("%8.2f", value);   % -1000 >x>-10000
        elseif value > -100000.
            field = sprintf("%8.1f", value);   % -10000>x>-100000
        elseif value <= -999999.5
            field = print_scientific_8(value);
            return % field
        else
            field = sprintf("%8.1f", value);
            %try:
            % ifield = field.index('.');
            %except ValueError:
            %    raise ValueError('error printing float; cant find decimal; field=%r value=%s' % (
            %        field, value))
            ifields = strfind(field, '.');
            if length(ifields) == 0
                error('error printing float; cant find decimal; field=%r value=%s', field, num2str(value))
            end
            if ifield < 8
                ival = int(round(value, 0));
                field = sprintf('%7d.', ival);
                %assert '.' != field[0], field
            else
                field = print_scientific_8(value);
            end
            return % field
        end
    end
    % original = '  000 Hello World 000  ';
    % cleaned = regexprep(original, '^[ 0]+|[ 0]+$', '');
    % disp(cleaned);  % This will output: 'Hello World'

    %field = field.strip(' 0');  % python
    %field = strip(field, 'both', ' 0');  % 2 calls
    field = regexprep(field, '^[ 0]+|[ 0]+$', '');

    if not(isstring(field))
        error('not a string')
    end
    field = sprintf('%8s', field);

    %assert length(field) == 8, ('value=%r field=%r is not 8 characters '
    %                           'long, its %s' % (value, field, length(field)))
    return  % field
end

function field = print_scientific_8(value)
    % Prints a value in 8-character scientific notation.
    % This is a sub-method and shouldn't typically be called
    %
    % Notes
    % -----
    % print_float_8 : a better float printing method

    if value == 0.0
        field = sprintf('%8s', '0.');
        return  % field
    end

    python_value = sprintf('%8.11e', value);

    %svalue, sexponent = strip(python_value).split('e')
    svalue_sexponent = split(strip(python_value), 'e');
    svalue = svalue_sexponent{1};
    sexponent = svalue_sexponent{2};
    
    % exponent = int(sexponent);  % removes 0s
    exponent = str2double(sexponent);  % removes 0s

    if abs(value) < 1.
        sign = '-';
    else
        sign = '+';
    end
    %sign = '-' if abs(value) < 1. else '+'

    % the exponent will be added later...
    % exp2 = str(exponent).strip('-+');
    exp2 = regexprep(string(exponent), '^[-+]+|[-+]+$', '');
    
    %value2 = float(svalue);
    value2 = str2double(svalue);

    leftover = 5 - length(exp2);
    if value < 0
        fmt = sprintf("%%1.%df", leftover - 1);
    else
        fmt = sprintf("%%1.%df", leftover);
    end
    svalue3 = sprintf(fmt, value2);
    svalue4 = strip(svalue3, '0');
    field = sprintf("%8s", svalue4 + sign + exp2);
    return  % field
end
