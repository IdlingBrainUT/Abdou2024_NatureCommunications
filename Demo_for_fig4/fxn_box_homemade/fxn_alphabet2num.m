%% conversion alphabet input to number
function convert_num = fxn_alphabet2num(input_str)
%% for debug,
% input_str = ('B');
%%
if input_str == ('A')
convert_num = 1;
elseif input_str == ('B')
convert_num = 2;
elseif input_str == ('C')
convert_num = 3;
elseif input_str == ('D')
convert_num = 4;
elseif input_str == ('E')
convert_num = 5;
elseif input_str == ('F')
convert_num = 6;
elseif input_str == ('G')
convert_num = 7;
elseif input_str == ('H')
convert_num = 8;
elseif input_str == ('I')
convert_num = 9;
elseif input_str == ('J')
convert_num = 10;
elseif input_str == ('K')
convert_num = 11;
elseif input_str == ('L')
convert_num = 12;
elseif input_str == ('M')
convert_num = 13;
elseif input_str == ('N')
convert_num = 14;
elseif input_str == ('O')
convert_num = 15;
elseif input_str == ('P')
convert_num = 16;
elseif input_str == ('Q')
convert_num = 17;
elseif input_str == ('R')
convert_num = 18;
elseif input_str == ('S')
convert_num = 19;
elseif input_str == ('T')
convert_num = 20;
elseif input_str == ('U')
convert_num = 21;
elseif input_str == ('V')
convert_num = 22;
elseif input_str == ('W')
convert_num = 23;
elseif input_str == ('X')
convert_num = 24;
elseif input_str == ('Y')
convert_num = 25;
elseif input_str == ('Z')
convert_num = 26;
end
%%
% disp([input_str,(' is converted to '), num2str(convert_num)]);
%%
end