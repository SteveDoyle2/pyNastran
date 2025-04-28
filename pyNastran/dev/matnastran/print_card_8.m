function out = print_card_8(fields)
    % Prints a nastran-style card with 8-character width fields.
    %
    % Parameters
    % ----------
    % fields : list[int/float/str/None] / cell-array
    %     all the fields in the BDF card (no trailing Nones)
    %
    % Returns
    % -------
    % card : str
    %     string representation of the card in small field format
    %
    %.. note:: An internal field value of None or '' will be treated as
    %          a blank field
    %.. note:: A small field format follows the  8-8-8-8-8-8-8-8 = 80
    %          format where the first 8 is the card name or
    %          blank (continuation).  The last 8-character field indicates
    %          an optional continuation, but because it's a left-justified
    %          unnecessary field, print_card doesn't use it.
    %
    %.. code-block:: python
    %
    %   >>> fields = ['DUMMY', 1, 2, 3, None, 4, 5, 6, 7, 8.]
    %   >>> fields = {'DUMMY', 1, 2, 3, None, 4, 5, 6, 7, 8.}
    %   >>> print_card_8(fields)
    %   DUMMY          1       2       3               4       5       6       7
    %                 8.
    %
    % print_card_8({'GRID', int32(1), int32(2), 3.0})
    fields
    %try
        name = fields(1);
        out = sprintf('%-8s', name{1})
    %catch
    %    warning("ERROR!  fields=%s"); % % fields)
        %sys.stdout.flush()
        %raise
    %end

    for i=2:length(fields)
        fielda = fields(i);
        field = fielda{1};
        disp(sprintf('field[%d] = %s', i, num2str(field)));
        %try
            %out += print_field_8(field);
            fieldi = print_field_8(field)
            out = strcat(out, fieldi)
        %catch
        %    fields
        %    warning("bad fields");
        %    %raise
        %end
        if mod(i-1, 8) == 0  % allow 1+8 fields per line
            %out = out.rstrip(' ')
            out
            outa = strip(out, 'right');
            %if out[-1] == '\n'  # empty line
            outi = outa(end);
            if strcmp(outa(end), '\n')  % empty line
                out = strcat(out, '+')
            end
            out = strcat(out, '\n        ');
        end
    end
    % removes blank lines at the end of cards
    out = regexprep(out, '[ \n+]+$', '');
    % out = out.rstrip(' \n+') + '\n'  
    return % out
end

function field = print_field_8(value)
    % Prints an 8-character width field
    %
    % Parameters
    % ----------
    % value : int/float/str
    %     the value to print
    % 
    % Returns
    % -------
    % field : str
    %     an 8-character string

    if isinteger(value)
        field = sprintf('%8d', value);
    elseif isnumeric(value)
        field = print_float_8(value);
    %elseif value is None
        %field = '        ';
    else
        field = sprintf('%8s', value);
    end
    if length(field) ~= 8
        error("field='%s' is not 8 characters long...raw_value=%s", field, num2str(value));
    end
    return  % field
end
