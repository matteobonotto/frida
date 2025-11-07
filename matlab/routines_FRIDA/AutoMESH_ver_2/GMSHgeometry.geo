ii=1

ii =

     1

element=shape{ii};
    spacing=refinement{ii};
    for i=1:size(element,1)
        disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
            ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
        n_point=n_point+1;
    end
{Undefined function or variable 'n_point'.
} 
n_point=1;
element=shape{ii};
    spacing=refinement{ii};
    for i=1:size(element,1)
        disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
            ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
        n_point=n_point+1;
    end
Point(1)={3,0,0,0.1};
{Index exceeds matrix dimensions.
} 
disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
            ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
{Index exceeds matrix dimensions.
} 
num2str(element(i,1))

ans =

    '2.998'

num2str(element(i,2))

ans =

    '0.063424'

num2str(spacing(i))
{Index exceeds matrix dimensions.
} 
spacing(size(elements,1),1)=refinement{ii};
{Undefined function or variable 'elements'.
} 
element=shape{ii};
spacing(size(element,1),1)=refinement{ii};
spacing=zeros(size(element,1),1);
spacing(:)=refinement{ii};
for i=1:size(element,1)
        disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
            ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
        n_point=n_point+1;
    end
Point(2)={3,0,0,0.1};
Point(3)={2.998,0.063424,0,0.1};
Point(4)={2.992,0.12659,0,0.1};
Point(5)={2.9819,0.18925,0,0.1};
Point(6)={2.9679,0.25115,0,0.1};
Point(7)={2.9501,0.31203,0,0.1};
Point(8)={2.9284,0.37166,0,0.1};
Point(9)={2.9029,0.42979,0,0.1};
Point(10)={2.8738,0.4862,0,0.1};
Point(11)={2.8413,0.54064,0,0.1};
Point(12)={2.8053,0.59291,0,0.1};
Point(13)={2.766,0.64279,0,0.1};
Point(14)={2.7237,0.69008,0,0.1};
Point(15)={2.6785,0.73459,0,0.1};
Point(16)={2.6306,0.77615,0,0.1};
Point(17)={2.5801,0.81458,0,0.1};
Point(18)={2.5272,0.84973,0,0.1};
Point(19)={2.4723,0.88145,0,0.1};
Point(20)={2.4154,0.90963,0,0.1};
Point(21)={2.3569,0.93415,0,0.1};
Point(22)={2.2969,0.9549,0,0.1};
Point(23)={2.2358,0.97181,0,0.1};
Point(24)={2.1736,0.98481,0,0.1};
Point(25)={2.1108,0.99384,0,0.1};
Point(26)={2.0476,0.99887,0,0.1};
Point(27)={1.9841,0.99987,0,0.1};
Point(28)={1.9208,0.99685,0,0.1};
Point(29)={1.8577,0.98982,0,0.1};
Point(30)={1.7952,0.9788,0,0.1};
Point(31)={1.7335,0.96384,0,0.1};
Point(32)={1.6729,0.945,0,0.1};
Point(33)={1.6137,0.92235,0,0.1};
Point(34)={1.5559,0.89599,0,0.1};
Point(35)={1.5,0.86603,0,0.1};
Point(36)={1.4461,0.83257,0,0.1};
Point(37)={1.3944,0.79576,0,0.1};
Point(38)={1.3451,0.75575,0,0.1};
Point(39)={1.2985,0.71269,0,0.1};
Point(40)={1.2547,0.66677,0,0.1};
Point(41)={1.2139,0.61816,0,0.1};
Point(42)={1.1763,0.56706,0,0.1};
Point(43)={1.142,0.51368,0,0.1};
Point(44)={1.1112,0.45823,0,0.1};
Point(45)={1.0839,0.40093,0,0.1};
Point(46)={1.0603,0.34202,0,0.1};
Point(47)={1.0405,0.28173,0,0.1};
Point(48)={1.0246,0.22031,0,0.1};
Point(49)={1.0126,0.158,0,0.1};
Point(50)={1.0045,0.095056,0,0.1};
Point(51)={1.0005,0.031728,0,0.1};
Point(52)={1.0005,-0.031728,0,0.1};
Point(53)={1.0045,-0.095056,0,0.1};
Point(54)={1.0126,-0.158,0,0.1};
Point(55)={1.0246,-0.22031,0,0.1};
Point(56)={1.0405,-0.28173,0,0.1};
Point(57)={1.0603,-0.34202,0,0.1};
Point(58)={1.0839,-0.40093,0,0.1};
Point(59)={1.1112,-0.45823,0,0.1};
Point(60)={1.142,-0.51368,0,0.1};
Point(61)={1.1763,-0.56706,0,0.1};
Point(62)={1.2139,-0.61816,0,0.1};
Point(63)={1.2547,-0.66677,0,0.1};
Point(64)={1.2985,-0.71269,0,0.1};
Point(65)={1.3451,-0.75575,0,0.1};
Point(66)={1.3944,-0.79576,0,0.1};
Point(67)={1.4461,-0.83257,0,0.1};
Point(68)={1.5,-0.86603,0,0.1};
Point(69)={1.5559,-0.89599,0,0.1};
Point(70)={1.6137,-0.92235,0,0.1};
Point(71)={1.6729,-0.945,0,0.1};
Point(72)={1.7335,-0.96384,0,0.1};
Point(73)={1.7952,-0.9788,0,0.1};
Point(74)={1.8577,-0.98982,0,0.1};
Point(75)={1.9208,-0.99685,0,0.1};
Point(76)={1.9841,-0.99987,0,0.1};
Point(77)={2.0476,-0.99887,0,0.1};
Point(78)={2.1108,-0.99384,0,0.1};
Point(79)={2.1736,-0.98481,0,0.1};
Point(80)={2.2358,-0.97181,0,0.1};
Point(81)={2.2969,-0.9549,0,0.1};
Point(82)={2.3569,-0.93415,0,0.1};
Point(83)={2.4154,-0.90963,0,0.1};
Point(84)={2.4723,-0.88145,0,0.1};
Point(85)={2.5272,-0.84973,0,0.1};
Point(86)={2.5801,-0.81458,0,0.1};
Point(87)={2.6306,-0.77615,0,0.1};
Point(88)={2.6785,-0.73459,0,0.1};
Point(89)={2.7237,-0.69008,0,0.1};
Point(90)={2.766,-0.64279,0,0.1};
Point(91)={2.8053,-0.59291,0,0.1};
Point(92)={2.8413,-0.54064,0,0.1};
Point(93)={2.8738,-0.4862,0,0.1};
Point(94)={2.9029,-0.42979,0,0.1};
Point(95)={2.9284,-0.37166,0,0.1};
Point(96)={2.9501,-0.31203,0,0.1};
Point(97)={2.9679,-0.25115,0,0.1};
Point(98)={2.9819,-0.18925,0,0.1};
Point(99)={2.992,-0.12659,0,0.1};
Point(100)={2.998,-0.063424,0,0.1};
Point(101)={3,-2.4493e-16,0,0.1};
for i=1:size(element,1)
        if i<size(element,1)
            disp(['Line(',num2str(n_line),')={', num2str(n_line),...
                ',', num2str(n_line+1), '};'])
        elseif i==size(element,1)
            disp(['Line(',num2str(n_line),')={', num2str(n_line),...
                ',', num2str(n_line-size(element,1)+1), '};'])
        end
        n_line=n_line+1;
    end
{Undefined function or variable 'n_line'.
} 
n_point=1;
n_line=1;
n_loop=1;
element=shape{ii};
    spacing=zeros(size(element,1),1);
    spacing(:)=refinement{ii};
    for i=1:size(element,1)
        disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
            ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
        n_point=n_point+1;
    end
Point(1)={3,0,0,0.1};
Point(2)={2.998,0.063424,0,0.1};
Point(3)={2.992,0.12659,0,0.1};
Point(4)={2.9819,0.18925,0,0.1};
Point(5)={2.9679,0.25115,0,0.1};
Point(6)={2.9501,0.31203,0,0.1};
Point(7)={2.9284,0.37166,0,0.1};
Point(8)={2.9029,0.42979,0,0.1};
Point(9)={2.8738,0.4862,0,0.1};
Point(10)={2.8413,0.54064,0,0.1};
Point(11)={2.8053,0.59291,0,0.1};
Point(12)={2.766,0.64279,0,0.1};
Point(13)={2.7237,0.69008,0,0.1};
Point(14)={2.6785,0.73459,0,0.1};
Point(15)={2.6306,0.77615,0,0.1};
Point(16)={2.5801,0.81458,0,0.1};
Point(17)={2.5272,0.84973,0,0.1};
Point(18)={2.4723,0.88145,0,0.1};
Point(19)={2.4154,0.90963,0,0.1};
Point(20)={2.3569,0.93415,0,0.1};
Point(21)={2.2969,0.9549,0,0.1};
Point(22)={2.2358,0.97181,0,0.1};
Point(23)={2.1736,0.98481,0,0.1};
Point(24)={2.1108,0.99384,0,0.1};
Point(25)={2.0476,0.99887,0,0.1};
Point(26)={1.9841,0.99987,0,0.1};
Point(27)={1.9208,0.99685,0,0.1};
Point(28)={1.8577,0.98982,0,0.1};
Point(29)={1.7952,0.9788,0,0.1};
Point(30)={1.7335,0.96384,0,0.1};
Point(31)={1.6729,0.945,0,0.1};
Point(32)={1.6137,0.92235,0,0.1};
Point(33)={1.5559,0.89599,0,0.1};
Point(34)={1.5,0.86603,0,0.1};
Point(35)={1.4461,0.83257,0,0.1};
Point(36)={1.3944,0.79576,0,0.1};
Point(37)={1.3451,0.75575,0,0.1};
Point(38)={1.2985,0.71269,0,0.1};
Point(39)={1.2547,0.66677,0,0.1};
Point(40)={1.2139,0.61816,0,0.1};
Point(41)={1.1763,0.56706,0,0.1};
Point(42)={1.142,0.51368,0,0.1};
Point(43)={1.1112,0.45823,0,0.1};
Point(44)={1.0839,0.40093,0,0.1};
Point(45)={1.0603,0.34202,0,0.1};
Point(46)={1.0405,0.28173,0,0.1};
Point(47)={1.0246,0.22031,0,0.1};
Point(48)={1.0126,0.158,0,0.1};
Point(49)={1.0045,0.095056,0,0.1};
Point(50)={1.0005,0.031728,0,0.1};
Point(51)={1.0005,-0.031728,0,0.1};
Point(52)={1.0045,-0.095056,0,0.1};
Point(53)={1.0126,-0.158,0,0.1};
Point(54)={1.0246,-0.22031,0,0.1};
Point(55)={1.0405,-0.28173,0,0.1};
Point(56)={1.0603,-0.34202,0,0.1};
Point(57)={1.0839,-0.40093,0,0.1};
Point(58)={1.1112,-0.45823,0,0.1};
Point(59)={1.142,-0.51368,0,0.1};
Point(60)={1.1763,-0.56706,0,0.1};
Point(61)={1.2139,-0.61816,0,0.1};
Point(62)={1.2547,-0.66677,0,0.1};
Point(63)={1.2985,-0.71269,0,0.1};
Point(64)={1.3451,-0.75575,0,0.1};
Point(65)={1.3944,-0.79576,0,0.1};
Point(66)={1.4461,-0.83257,0,0.1};
Point(67)={1.5,-0.86603,0,0.1};
Point(68)={1.5559,-0.89599,0,0.1};
Point(69)={1.6137,-0.92235,0,0.1};
Point(70)={1.6729,-0.945,0,0.1};
Point(71)={1.7335,-0.96384,0,0.1};
Point(72)={1.7952,-0.9788,0,0.1};
Point(73)={1.8577,-0.98982,0,0.1};
Point(74)={1.9208,-0.99685,0,0.1};
Point(75)={1.9841,-0.99987,0,0.1};
Point(76)={2.0476,-0.99887,0,0.1};
Point(77)={2.1108,-0.99384,0,0.1};
Point(78)={2.1736,-0.98481,0,0.1};
Point(79)={2.2358,-0.97181,0,0.1};
Point(80)={2.2969,-0.9549,0,0.1};
Point(81)={2.3569,-0.93415,0,0.1};
Point(82)={2.4154,-0.90963,0,0.1};
Point(83)={2.4723,-0.88145,0,0.1};
Point(84)={2.5272,-0.84973,0,0.1};
Point(85)={2.5801,-0.81458,0,0.1};
Point(86)={2.6306,-0.77615,0,0.1};
Point(87)={2.6785,-0.73459,0,0.1};
Point(88)={2.7237,-0.69008,0,0.1};
Point(89)={2.766,-0.64279,0,0.1};
Point(90)={2.8053,-0.59291,0,0.1};
Point(91)={2.8413,-0.54064,0,0.1};
Point(92)={2.8738,-0.4862,0,0.1};
Point(93)={2.9029,-0.42979,0,0.1};
Point(94)={2.9284,-0.37166,0,0.1};
Point(95)={2.9501,-0.31203,0,0.1};
Point(96)={2.9679,-0.25115,0,0.1};
Point(97)={2.9819,-0.18925,0,0.1};
Point(98)={2.992,-0.12659,0,0.1};
Point(99)={2.998,-0.063424,0,0.1};
Point(100)={3,-2.4493e-16,0,0.1};
    
    disp(' ')
 
for i=1:size(element,1)
        if i<size(element,1)
            disp(['Line(',num2str(n_line),')={', num2str(n_line),...
                ',', num2str(n_line+1), '};'])
        elseif i==size(element,1)
            disp(['Line(',num2str(n_line),')={', num2str(n_line),...
                ',', num2str(n_line-size(element,1)+1), '};'])
        end
        n_line=n_line+1;
    end
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,9};
Line(9)={9,10};
Line(10)={10,11};
Line(11)={11,12};
Line(12)={12,13};
Line(13)={13,14};
Line(14)={14,15};
Line(15)={15,16};
Line(16)={16,17};
Line(17)={17,18};
Line(18)={18,19};
Line(19)={19,20};
Line(20)={20,21};
Line(21)={21,22};
Line(22)={22,23};
Line(23)={23,24};
Line(24)={24,25};
Line(25)={25,26};
Line(26)={26,27};
Line(27)={27,28};
Line(28)={28,29};
Line(29)={29,30};
Line(30)={30,31};
Line(31)={31,32};
Line(32)={32,33};
Line(33)={33,34};
Line(34)={34,35};
Line(35)={35,36};
Line(36)={36,37};
Line(37)={37,38};
Line(38)={38,39};
Line(39)={39,40};
Line(40)={40,41};
Line(41)={41,42};
Line(42)={42,43};
Line(43)={43,44};
Line(44)={44,45};
Line(45)={45,46};
Line(46)={46,47};
Line(47)={47,48};
Line(48)={48,49};
Line(49)={49,50};
Line(50)={50,51};
Line(51)={51,52};
Line(52)={52,53};
Line(53)={53,54};
Line(54)={54,55};
Line(55)={55,56};
Line(56)={56,57};
Line(57)={57,58};
Line(58)={58,59};
Line(59)={59,60};
Line(60)={60,61};
Line(61)={61,62};
Line(62)={62,63};
Line(63)={63,64};
Line(64)={64,65};
Line(65)={65,66};
Line(66)={66,67};
Line(67)={67,68};
Line(68)={68,69};
Line(69)={69,70};
Line(70)={70,71};
Line(71)={71,72};
Line(72)={72,73};
Line(73)={73,74};
Line(74)={74,75};
Line(75)={75,76};
Line(76)={76,77};
Line(77)={77,78};
Line(78)={78,79};
Line(79)={79,80};
Line(80)={80,81};
Line(81)={81,82};
Line(82)={82,83};
Line(83)={83,84};
Line(84)={84,85};
Line(85)={85,86};
Line(86)={86,87};
Line(87)={87,88};
Line(88)={88,89};
Line(89)={89,90};
Line(90)={90,91};
Line(91)={91,92};
Line(92)={92,93};
Line(93)={93,94};
Line(94)={94,95};
Line(95)={95,96};
Line(96)={96,97};
Line(97)={97,98};
Line(98)={98,99};
Line(99)={99,100};
Line(100)={100,1};
disp(' ')
 
    
    % Line Loop
    disp(['Line Loop(',num2str(n_loop),')={'])
Line Loop(1)={
    for i=1:size(element,1)
        if i<size(element,1)
            disp([num2str(n_line-size(element,1)+i-1), ','])
        elseif i==size(element,1)
            disp(num2str(n_line-size(element,1)+i-1))
        end
    end
1,
2,
3,
4,
5,
6,
7,
8,
9,
10,
11,
12,
13,
14,
15,
16,
17,
18,
19,
20,
21,
22,
23,
24,
25,
26,
27,
28,
29,
30,
31,
32,
33,
34,
35,
36,
37,
38,
39,
40,
41,
42,
43,
44,
45,
46,
47,
48,
49,
50,
51,
52,
53,
54,
55,
56,
57,
58,
59,
60,
61,
62,
63,
64,
65,
66,
67,
68,
69,
70,
71,
72,
73,
74,
75,
76,
77,
78,
79,
80,
81,
82,
83,
84,
85,
86,
87,
88,
89,
90,
91,
92,
93,
94,
95,
96,
97,
98,
99,
100
    n_loop=n_loop+1;
    disp('};')
};
    
    disp(' ')
 
element=shape{ii};
    spacing=zeros(size(element,1),1);
    spacing(:)=refinement{ii};
    for i=1:size(element,1)
        disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
            ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
        fprintf(fileID,'%6s 12.12f\n %12s\n','Point(', num2str(n_point) ,'exp(x)');
% %         fprintf(fileID,'%12.12f\n', ',' ,RR, '%6.2f\n', ',' , ZZ);

        n_point=n_point+1;
    end
Point(101)={3,0,0,0.1};
{Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('fprintf')" style="font-weight:bold">fprintf</a>
Invalid file identifier.  Use fopen to generate a
valid file identifier.
} 
fileID = fopen('GMSH_GEOMETRY.geo','w');
fprintf(fileID,'%12.12f\n', ',' ,RR, '%6.2f\n', ',' , ZZ);
fclose(fileID);
element=shape{ii};
    spacing=zeros(size(element,1),1);
    spacing(:)=refinement{ii};
    for i=1:size(element,1)
        disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
            ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
        fprintf(fileID,'%6s 12.12f\n %12s\n','Point(', num2str(n_point) ,'exp(x)');
% %         fprintf(fileID,'%12.12f\n', ',' ,RR, '%6.2f\n', ',' , ZZ);

        n_point=n_point+1;
    end
Point(101)={3,0,0,0.1};
{Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('fprintf')" style="font-weight:bold">fprintf</a>
Invalid file identifier.  Use fopen to generate a
valid file identifier.
} 
fileID

fileID =

    11

if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
{Undefined function or variable 'ii'.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('main_MESH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m', 21)" style="font-weight:bold">main_MESH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m',21,0)">line 21</a>)
element=shape{ii};
} 
opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m',21,0)
main_MESH
for i=1:size(element,1)
    disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
        ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
    fprintf(fileID,'%6s 12.12f\n %12s\n','Point(', num2str(n_point) ,'exp(x)');
    % %         fprintf(fileID,'%12.12f\n', ',' ,RR, '%6.2f\n', ',' , ZZ);
    
    n_point=n_point+1;
end
Point(1)={3,0,0,0.1};
{Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('fprintf')" style="font-weight:bold">fprintf</a>
Invalid file identifier.  Use fopen to generate a
valid file identifier.
} 
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
for i=1:size(element,1)
    disp(['Point(',num2str(n_point),')={', num2str(element(i,1)),...
        ',', num2str(element(i,2)), ',0,' ,num2str(spacing(i)), '};'])
    fprintf(fileID,'%6s 12.12f\n %12s\n','Point(', num2str(n_point) ,'exp(x)');
    % %         fprintf(fileID,'%12.12f\n', ',' ,RR, '%6.2f\n', ',' , ZZ);
    
    n_point=n_point+1;
end
Point(1)={3,0,0,0.1};
Point(2)={2.998,0.063424,0,0.1};
Point(3)={2.992,0.12659,0,0.1};
Point(4)={2.9819,0.18925,0,0.1};
Point(5)={2.9679,0.25115,0,0.1};
Point(6)={2.9501,0.31203,0,0.1};
Point(7)={2.9284,0.37166,0,0.1};
Point(8)={2.9029,0.42979,0,0.1};
Point(9)={2.8738,0.4862,0,0.1};
Point(10)={2.8413,0.54064,0,0.1};
Point(11)={2.8053,0.59291,0,0.1};
Point(12)={2.766,0.64279,0,0.1};
Point(13)={2.7237,0.69008,0,0.1};
Point(14)={2.6785,0.73459,0,0.1};
Point(15)={2.6306,0.77615,0,0.1};
Point(16)={2.5801,0.81458,0,0.1};
Point(17)={2.5272,0.84973,0,0.1};
Point(18)={2.4723,0.88145,0,0.1};
Point(19)={2.4154,0.90963,0,0.1};
Point(20)={2.3569,0.93415,0,0.1};
Point(21)={2.2969,0.9549,0,0.1};
Point(22)={2.2358,0.97181,0,0.1};
Point(23)={2.1736,0.98481,0,0.1};
Point(24)={2.1108,0.99384,0,0.1};
Point(25)={2.0476,0.99887,0,0.1};
Point(26)={1.9841,0.99987,0,0.1};
Point(27)={1.9208,0.99685,0,0.1};
Point(28)={1.8577,0.98982,0,0.1};
Point(29)={1.7952,0.9788,0,0.1};
Point(30)={1.7335,0.96384,0,0.1};
Point(31)={1.6729,0.945,0,0.1};
Point(32)={1.6137,0.92235,0,0.1};
Point(33)={1.5559,0.89599,0,0.1};
Point(34)={1.5,0.86603,0,0.1};
Point(35)={1.4461,0.83257,0,0.1};
Point(36)={1.3944,0.79576,0,0.1};
Point(37)={1.3451,0.75575,0,0.1};
Point(38)={1.2985,0.71269,0,0.1};
Point(39)={1.2547,0.66677,0,0.1};
Point(40)={1.2139,0.61816,0,0.1};
Point(41)={1.1763,0.56706,0,0.1};
Point(42)={1.142,0.51368,0,0.1};
Point(43)={1.1112,0.45823,0,0.1};
Point(44)={1.0839,0.40093,0,0.1};
Point(45)={1.0603,0.34202,0,0.1};
Point(46)={1.0405,0.28173,0,0.1};
Point(47)={1.0246,0.22031,0,0.1};
Point(48)={1.0126,0.158,0,0.1};
Point(49)={1.0045,0.095056,0,0.1};
Point(50)={1.0005,0.031728,0,0.1};
Point(51)={1.0005,-0.031728,0,0.1};
Point(52)={1.0045,-0.095056,0,0.1};
Point(53)={1.0126,-0.158,0,0.1};
Point(54)={1.0246,-0.22031,0,0.1};
Point(55)={1.0405,-0.28173,0,0.1};
Point(56)={1.0603,-0.34202,0,0.1};
Point(57)={1.0839,-0.40093,0,0.1};
Point(58)={1.1112,-0.45823,0,0.1};
Point(59)={1.142,-0.51368,0,0.1};
Point(60)={1.1763,-0.56706,0,0.1};
Point(61)={1.2139,-0.61816,0,0.1};
Point(62)={1.2547,-0.66677,0,0.1};
Point(63)={1.2985,-0.71269,0,0.1};
Point(64)={1.3451,-0.75575,0,0.1};
Point(65)={1.3944,-0.79576,0,0.1};
Point(66)={1.4461,-0.83257,0,0.1};
Point(67)={1.5,-0.86603,0,0.1};
Point(68)={1.5559,-0.89599,0,0.1};
Point(69)={1.6137,-0.92235,0,0.1};
Point(70)={1.6729,-0.945,0,0.1};
Point(71)={1.7335,-0.96384,0,0.1};
Point(72)={1.7952,-0.9788,0,0.1};
Point(73)={1.8577,-0.98982,0,0.1};
Point(74)={1.9208,-0.99685,0,0.1};
Point(75)={1.9841,-0.99987,0,0.1};
Point(76)={2.0476,-0.99887,0,0.1};
Point(77)={2.1108,-0.99384,0,0.1};
Point(78)={2.1736,-0.98481,0,0.1};
Point(79)={2.2358,-0.97181,0,0.1};
Point(80)={2.2969,-0.9549,0,0.1};
Point(81)={2.3569,-0.93415,0,0.1};
Point(82)={2.4154,-0.90963,0,0.1};
Point(83)={2.4723,-0.88145,0,0.1};
Point(84)={2.5272,-0.84973,0,0.1};
Point(85)={2.5801,-0.81458,0,0.1};
Point(86)={2.6306,-0.77615,0,0.1};
Point(87)={2.6785,-0.73459,0,0.1};
Point(88)={2.7237,-0.69008,0,0.1};
Point(89)={2.766,-0.64279,0,0.1};
Point(90)={2.8053,-0.59291,0,0.1};
Point(91)={2.8413,-0.54064,0,0.1};
Point(92)={2.8738,-0.4862,0,0.1};
Point(93)={2.9029,-0.42979,0,0.1};
Point(94)={2.9284,-0.37166,0,0.1};
Point(95)={2.9501,-0.31203,0,0.1};
Point(96)={2.9679,-0.25115,0,0.1};
Point(97)={2.9819,-0.18925,0,0.1};
Point(98)={2.992,-0.12659,0,0.1};
Point(99)={2.998,-0.063424,0,0.1};
Point(100)={3,-2.4493e-16,0,0.1};
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
P=[6.5443  871.2282];
c=[6.7469 795.1066];
c1=[6.5472 870.5470];
c2=[6.4626 910.2019];
table=[P;c;c1;c2];
words={'polyfit';'basis functions';'nlinfit';'L1'};
stuff={table,words};
for k1 = 1:size(words,1)
    fprintf('%.4f \t %.4f \t %s \n',stuff{1}(k1,:),char(stuff{2}(k1)))
end
6.5443 	 871.2282 	 polyfit 
6.7469 	 795.1066 	 basis functions 
6.5472 	 870.5470 	 nlinfit 
6.4626 	 910.2019 	 L1 
P=[6.5443  871.2282];
c=[6.7469 795.1066];
c1=[6.5472 870.5470];
c2=[6.4626 910.2019];
table=[P;c;c1;c2];
words={'polyfit';'basis functions';'nlinfit';'L1'};
stuff={table,words};
for k1 = 1:size(words,1)
    fprintf(fileID,'%.4f \t %.4f \t %s \n',stuff{1}(k1,:),char(stuff{2}(k1)))
end

ans =

    29


ans =

    37


ans =

    29


ans =

    24

char(stuff{2}(k1))

ans =

    'L1'

if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
for k1 = 1:size(words,1)
    fprintf(fileID,'%s  \t %.4f \t %.4f \t %s \n',char('poly'),stuff{1}(k1,:),char(stuff{2}(k1)))
end

ans =

    37


ans =

    45


ans =

    37


ans =

    32

fclose(fileID);
stuff{1}

ans =

   1.0e+02 *

   0.065443000000000   8.712282000000000
   0.067469000000000   7.951066000000000
   0.065472000000000   8.705470000000000
   0.064626000000000   9.102019000000000

if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH

ans =

    38


ans =

    47

{Index exceeds matrix dimensions.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('main_MESH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m', 45)" style="font-weight:bold">main_MESH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m',45,0)">line 45</a>)
    fprintf(fileID,'%s  \t %.4f \t %.4f \t %s
    \n',char('Point('),c1(k1),...
} 
opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m',45,0)
c1

c1 =

   1.0e+02 *

   0.065472000000000   8.705470000000000

main_MESH

ans =

    38


ans =

    45


ans =

    38


ans =

    32

if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH

ans =

    36


ans =

    43


ans =

    36


ans =

    30

if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH

ans =

    36


ans =

    43


ans =

    36


ans =

    30

if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
Exception in thread "AWT-EventQueue-0" java.lang.NullPointerException
	at org.netbeans.editor.BaseDocument.notifyUnmodify(BaseDocument.java:1465)
	at org.netbeans.editor.BaseDocument.notifyModifyCheckEnd(BaseDocument.java:816)
	at org.netbeans.editor.BaseDocumentEvent.redo(BaseDocumentEvent.java:336)
	at javax.swing.undo.UndoManager.redoTo(Unknown Source)
	at javax.swing.undo.UndoManager.redo(Unknown Source)
	at com.mathworks.mwswing.undo.MUndoManager.redo(MUndoManager.java:255)
	at org.netbeans.editor.ActionFactory$RedoAction.actionPerformed(ActionFactory.java:767)
	at org.netbeans.editor.BaseAction.actionPerformed(BaseAction.java:259)
	at javax.swing.SwingUtilities.notifyAction(Unknown Source)
	at javax.swing.JComponent.processKeyBinding(Unknown Source)
	at javax.swing.JComponent.processKeyBindings(Unknown Source)
	at javax.swing.JComponent.processKeyEvent(Unknown Source)
	at com.mathworks.widgets.SyntaxTextPaneBase.processKeyEvent(SyntaxTextPaneBase.java:1189)
	at java.awt.Component.processEvent(Unknown Source)
	at java.awt.Container.processEvent(Unknown Source)
	at java.awt.Component.dispatchEventImpl(Unknown Source)
	at java.awt.Container.dispatchEventImpl(Unknown Source)
	at java.awt.Component.dispatchEvent(Unknown Source)
	at java.awt.KeyboardFocusManager.redispatchEvent(Unknown Source)
	at java.awt.DefaultKeyboardFocusManager.dispatchKeyEvent(Unknown Source)
	at java.awt.DefaultKeyboardFocusManager.preDispatchKeyEvent(Unknown Source)
	at java.awt.DefaultKeyboardFocusManager.typeAheadAssertions(Unknown Source)
	at java.awt.DefaultKeyboardFocusManager.dispatchEvent(Unknown Source)
	at java.awt.Component.dispatchEventImpl(Unknown Source)
	at java.awt.Container.dispatchEventImpl(Unknown Source)
	at java.awt.Window.dispatchEventImpl(Unknown Source)
	at java.awt.Component.dispatchEvent(Unknown Source)
	at java.awt.EventQueue.dispatchEventImpl(Unknown Source)
	at java.awt.EventQueue.access$200(Unknown Source)
	at java.awt.EventQueue$3.run(Unknown Source)
	at java.awt.EventQueue$3.run(Unknown Source)
	at java.security.AccessController.doPrivileged(Native Method)
	at java.security.ProtectionDomain$1.doIntersectionPrivilege(Unknown Source)
	at java.security.ProtectionDomain$1.doIntersectionPrivilege(Unknown Source)
	at java.awt.EventQueue$4.run(Unknown Source)
	at java.awt.EventQueue$4.run(Unknown Source)
	at java.security.AccessController.doPrivileged(Native Method)
	at java.security.ProtectionDomain$1.doIntersectionPrivilege(Unknown Source)
	at java.awt.EventQueue.dispatchEvent(Unknown Source)
	at java.awt.EventDispatchThread.pumpOneEventForFilters(Unknown Source)
	at java.awt.EventDispatchThread.pumpEventsForFilter(Unknown Source)
	at java.awt.EventDispatchThread.pumpEventsForHierarchy(Unknown Source)
	at java.awt.EventDispatchThread.pumpEvents(Unknown Source)
	at java.awt.EventDispatchThread.pumpEvents(Unknown Source)
	at java.awt.EventDispatchThread.run(Unknown Source)
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
for ii = 1:size(element,1)
    fprintf(fileID,'%%d %s\n',...
        n_line-size(element,1)+i-1,...
        ',');
end
n_loop=n_loop+1;
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
{Undefined function or variable 'GMSH_GEOMETRY'.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('main_MESH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m', 14)" style="font-weight:bold">main_MESH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m',14,0)">line 14</a>)
filename=GMSH_GEOMETRY;
} 
main_MESH

fileID =

    34

[mesh]=readGMSH(filename);
{Index exceeds matrix dimensions.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('dlmread', 'C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m', 159)" style="font-weight:bold">dlmread</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m',159,0)">line 159</a>)
        result= result(:,1:ncols);

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('readGMSH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\readGMSH.m', 11)" style="font-weight:bold">readGMSH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\readGMSH.m',11,0)">line 11</a>)
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);
} 
main_MESH
file=([filename '.msh']);

% no of nodes is mentioned in 5th row and first column
N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);

node_id     = dlmread(file,'',[5 0 4+N_n 0]);
nodes       = dlmread(file,'',[5 1 4+N_n 3]);
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);
{Index exceeds matrix dimensions.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('dlmread', 'C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m', 159)" style="font-weight:bold">dlmread</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m',159,0)">line 159</a>)
        result= result(:,1:ncols);
} 
% no of nodes is mentioned in 5th row and first column
N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);
node_id     = dlmread(file,'',[5 0 4+N_n 0]);
nodes       = dlmread(file,'',[5 1 4+N_n 3]);
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);
{Index exceeds matrix dimensions.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('dlmread', 'C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m', 159)" style="font-weight:bold">dlmread</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m',159,0)">line 159</a>)
        result= result(:,1:ncols);
} 
uiopen('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\GMSH_GEOMETRY.geo',1)
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
Info    : Running 'gmsh -2 GMSH_GEOMETRY.geo' [Gmsh 2.16.0, 1 node, max. 1 thread]
Info    : Started on Mon Feb 26 15:37:36 2018
Info    : Reading 'GMSH_GEOMETRY.geo'...
Info    : Done reading 'GMSH_GEOMETRY.geo'
Info    : Finalized high order topology of periodic connections
Info    : Meshing 1D...
Info    : Meshing curve 1 (Line)
Info    : Meshing curve 2 (Line)
Info    : Meshing curve 3 (Line)
Info    : Meshing curve 4 (Line)
Info    : Meshing curve 5 (Line)
Info    : Meshing curve 6 (Line)
Info    : Meshing curve 7 (Line)
Info    : Meshing curve 8 (Line)
Info    : Meshing curve 9 (Line)
Info    : Meshing curve 10 (Line)
Info    : Meshing curve 11 (Line)
Info    : Meshing curve 12 (Line)
Info    : Meshing curve 13 (Line)
Info    : Meshing curve 14 (Line)
Info    : Meshing curve 15 (Line)
Info    : Meshing curve 16 (Line)
Info    : Meshing curve 17 (Line)
Info    : Meshing curve 18 (Line)
Info    : Meshing curve 19 (Line)
Info    : Meshing curve 20 (Line)
Info    : Meshing curve 21 (Line)
Info    : Meshing curve 22 (Line)
Info    : Meshing curve 23 (Line)
Info    : Meshing curve 24 (Line)
Info    : Meshing curve 25 (Line)
Info    : Meshing curve 26 (Line)
Info    : Meshing curve 27 (Line)
Info    : Meshing curve 28 (Line)
Info    : Meshing curve 29 (Line)
Info    : Meshing curve 30 (Line)
Info    : Meshing curve 31 (Line)
Info    : Meshing curve 32 (Line)
Info    : Meshing curve 33 (Line)
Info    : Meshing curve 34 (Line)
Info    : Meshing curve 35 (Line)
Info    : Meshing curve 36 (Line)
Info    : Meshing curve 37 (Line)
Info    : Meshing curve 38 (Line)
Info    : Meshing curve 39 (Line)
Info    : Meshing curve 40 (Line)
Info    : Meshing curve 41 (Line)
Info    : Meshing curve 42 (Line)
Info    : Meshing curve 43 (Line)
Info    : Meshing curve 44 (Line)
Info    : Meshing curve 45 (Line)
Info    : Meshing curve 46 (Line)
Info    : Meshing curve 47 (Line)
Info    : Meshing curve 48 (Line)
Info    : Meshing curve 49 (Line)
Info    : Meshing curve 50 (Line)
Info    : Meshing curve 51 (Line)
Info    : Meshing curve 52 (Line)
Info    : Meshing curve 53 (Line)
Info    : Meshing curve 54 (Line)
Info    : Meshing curve 55 (Line)
Info    : Meshing curve 56 (Line)
Info    : Meshing curve 57 (Line)
Info    : Meshing curve 58 (Line)
Info    : Meshing curve 59 (Line)
Info    : Meshing curve 60 (Line)
Info    : Meshing curve 61 (Line)
Info    : Meshing curve 62 (Line)
Info    : Meshing curve 63 (Line)
Info    : Meshing curve 64 (Line)
Info    : Meshing curve 65 (Line)
Info    : Meshing curve 66 (Line)
Info    : Meshing curve 67 (Line)
Info    : Meshing curve 68 (Line)
Info    : Meshing curve 69 (Line)
Info    : Meshing curve 70 (Line)
Info    : Meshing curve 71 (Line)
Info    : Meshing curve 72 (Line)
Info    : Meshing curve 73 (Line)
Info    : Meshing curve 74 (Line)
Info    : Meshing curve 75 (Line)
Info    : Meshing curve 76 (Line)
Info    : Meshing curve 77 (Line)
Info    : Meshing curve 78 (Line)
Info    : Meshing curve 79 (Line)
Info    : Meshing curve 80 (Line)
Info    : Meshing curve 81 (Line)
Info    : Meshing curve 82 (Line)
Info    : Meshing curve 83 (Line)
Info    : Meshing curve 84 (Line)
Info    : Meshing curve 85 (Line)
Info    : Meshing curve 86 (Line)
Info    : Meshing curve 87 (Line)
Info    : Meshing curve 88 (Line)
Info    : Meshing curve 89 (Line)
Info    : Meshing curve 90 (Line)
Info    : Meshing curve 91 (Line)
Info    : Meshing curve 92 (Line)
Info    : Meshing curve 93 (Line)
Info    : Meshing curve 94 (Line)
Info    : Meshing curve 95 (Line)
Info    : Meshing curve 96 (Line)
Info    : Meshing curve 97 (Line)
Info    : Meshing curve 98 (Line)
Info    : Meshing curve 99 (Line)
Info    : Meshing curve 100 (Line)
Warning : Curve 100 has a zero length
Info    : Done meshing 1D (0.015625 s)
Info    : Meshing 2D...
Info    : Done meshing 2D (0 s)
Info    : 100 vertices 200 elements
Warning : ------------------------------
Warning : Mesh generation error summary
Warning :     1 warning
Warning :     0 errors
Warning : Check the full log for details
Warning : ------------------------------
Info    : Writing 'GMSH_GEOMETRY.msh'...
Info    : Done writing 'GMSH_GEOMETRY.msh'
Info    : Stopped on Mon Feb 26 15:37:38 2018
{Index exceeds matrix dimensions.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('dlmread', 'C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m', 159)" style="font-weight:bold">dlmread</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m',159,0)">line 159</a>)
        result= result(:,1:ncols);

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('readGMSH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\readGMSH.m', 11)" style="font-weight:bold">readGMSH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\readGMSH.m',11,0)">line 11</a>)
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e
7]);

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('main_MESH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m', 18)" style="font-weight:bold">main_MESH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m',18,0)">line 18</a>)
[mesh]=readGMSH(filename);
} 
main_MESH
Info    : Running 'gmsh -2 GMSH_GEOMETRY.geo' [Gmsh 2.16.0, 1 node, max. 1 thread]
Info    : Started on Mon Feb 26 15:44:12 2018
Info    : Reading 'GMSH_GEOMETRY.geo'...
Error   : 'GMSH_GEOMETRY.geo', line 304 : syntax error ()
Info    : Done reading 'GMSH_GEOMETRY.geo'
Info    : Finalized high order topology of periodic connections
Info    : Meshing 1D...
Info    : Meshing curve 1 (Line)
Info    : Meshing curve 2 (Line)
Info    : Meshing curve 3 (Line)
Info    : Meshing curve 4 (Line)
Info    : Meshing curve 5 (Line)
Info    : Meshing curve 6 (Line)
Info    : Meshing curve 7 (Line)
Info    : Meshing curve 8 (Line)
Info    : Meshing curve 9 (Line)
Info    : Meshing curve 10 (Line)
Info    : Meshing curve 11 (Line)
Info    : Meshing curve 12 (Line)
Info    : Meshing curve 13 (Line)
Info    : Meshing curve 14 (Line)
Info    : Meshing curve 15 (Line)
Info    : Meshing curve 16 (Line)
Info    : Meshing curve 17 (Line)
Info    : Meshing curve 18 (Line)
Info    : Meshing curve 19 (Line)
Info    : Meshing curve 20 (Line)
Info    : Meshing curve 21 (Line)
Info    : Meshing curve 22 (Line)
Info    : Meshing curve 23 (Line)
Info    : Meshing curve 24 (Line)
Info    : Meshing curve 25 (Line)
Info    : Meshing curve 26 (Line)
Info    : Meshing curve 27 (Line)
Info    : Meshing curve 28 (Line)
Info    : Meshing curve 29 (Line)
Info    : Meshing curve 30 (Line)
Info    : Meshing curve 31 (Line)
Info    : Meshing curve 32 (Line)
Info    : Meshing curve 33 (Line)
Info    : Meshing curve 34 (Line)
Info    : Meshing curve 35 (Line)
Info    : Meshing curve 36 (Line)
Info    : Meshing curve 37 (Line)
Info    : Meshing curve 38 (Line)
Info    : Meshing curve 39 (Line)
Info    : Meshing curve 40 (Line)
Info    : Meshing curve 41 (Line)
Info    : Meshing curve 42 (Line)
Info    : Meshing curve 43 (Line)
Info    : Meshing curve 44 (Line)
Info    : Meshing curve 45 (Line)
Info    : Meshing curve 46 (Line)
Info    : Meshing curve 47 (Line)
Info    : Meshing curve 48 (Line)
Info    : Meshing curve 49 (Line)
Info    : Meshing curve 50 (Line)
Info    : Meshing curve 51 (Line)
Info    : Meshing curve 52 (Line)
Info    : Meshing curve 53 (Line)
Info    : Meshing curve 54 (Line)
Info    : Meshing curve 55 (Line)
Info    : Meshing curve 56 (Line)
Info    : Meshing curve 57 (Line)
Info    : Meshing curve 58 (Line)
Info    : Meshing curve 59 (Line)
Info    : Meshing curve 60 (Line)
Info    : Meshing curve 61 (Line)
Info    : Meshing curve 62 (Line)
Info    : Meshing curve 63 (Line)
Info    : Meshing curve 64 (Line)
Info    : Meshing curve 65 (Line)
Info    : Meshing curve 66 (Line)
Info    : Meshing curve 67 (Line)
Info    : Meshing curve 68 (Line)
Info    : Meshing curve 69 (Line)
Info    : Meshing curve 70 (Line)
Info    : Meshing curve 71 (Line)
Info    : Meshing curve 72 (Line)
Info    : Meshing curve 73 (Line)
Info    : Meshing curve 74 (Line)
Info    : Meshing curve 75 (Line)
Info    : Meshing curve 76 (Line)
Info    : Meshing curve 77 (Line)
Info    : Meshing curve 78 (Line)
Info    : Meshing curve 79 (Line)
Info    : Meshing curve 80 (Line)
Info    : Meshing curve 81 (Line)
Info    : Meshing curve 82 (Line)
Info    : Meshing curve 83 (Line)
Info    : Meshing curve 84 (Line)
Info    : Meshing curve 85 (Line)
Info    : Meshing curve 86 (Line)
Info    : Meshing curve 87 (Line)
Info    : Meshing curve 88 (Line)
Info    : Meshing curve 89 (Line)
Info    : Meshing curve 90 (Line)
Info    : Meshing curve 91 (Line)
Info    : Meshing curve 92 (Line)
Info    : Meshing curve 93 (Line)
Info    : Meshing curve 94 (Line)
Info    : Meshing curve 95 (Line)
Info    : Meshing curve 96 (Line)
Info    : Meshing curve 97 (Line)
Info    : Meshing curve 98 (Line)
Info    : Meshing curve 99 (Line)
Info    : Meshing curve 100 (Line)
Warning : Curve 100 has a zero length
Info    : Done meshing 1D (0 s)
Info    : Meshing 2D...
Info    : Done meshing 2D (0 s)
Info    : 100 vertices 200 elements
Warning : ------------------------------
Warning : Mesh generation error summary
Warning :     1 warning
Warning :     0 errors
Warning : Check the full log for details
Warning : ------------------------------
Info    : Writing 'GMSH_GEOMETRY.msh'...
Info    : Done writing 'GMSH_GEOMETRY.msh'
Info    : Stopped on Mon Feb 26 15:44:15 2018
{Index exceeds matrix dimensions.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('dlmread', 'C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m', 159)" style="font-weight:bold">dlmread</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\matlab\iofun\dlmread.m',159,0)">line 159</a>)
        result= result(:,1:ncols);

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('readGMSH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\readGMSH.m', 11)" style="font-weight:bold">readGMSH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\readGMSH.m',11,0)">line 11</a>)
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e
7]);

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('main_MESH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m', 18)" style="font-weight:bold">main_MESH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m',18,0)">line 18</a>)
[mesh]=readGMSH(filename);
} 
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
{Cell contents reference from a non-cell array object.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('write_GMSH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\write_GMSH.m', 12)" style="font-weight:bold">write_GMSH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\write_GMSH.m',12,0)">line 12</a>)
spacing(:)=refinement{1};

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('main_MESH', 'D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m', 12)" style="font-weight:bold">main_MESH</a> (<a href="matlab: opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\main_MESH.m',12,0)">line 12</a>)
[fileID] = write_GMSH(GEOMETRY,0.001,filename);
} 
opentoline('D:\Google Drive\PhD\RESEARCH_ACTIVITY\Mesh_Generator\gmsh-2.16.0\Sample_2\write_GMSH.m',12,0)
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
GEOMETRY([1 END],:)
{Undefined function or variable 'END'.
} 
GEOMETRY([1 end],:)

ans =

   3.000000000000000                   0
   3.000000000000000  -0.000000000000000

size(unique(GEOMETRY,2))
{Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('unique', 'C:\Program Files\MATLAB\R2017a\toolbox\matlab\ops\unique.p', 106)" style="font-weight:bold">unique</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\matlab\ops\unique.p',106,0)">line 106</a>)
Invalid input. Valid flags are 'rows', 'first',
'last', 'stable', 'sorted', 'legacy'.
} 
size(unique(GEOMETRY,rows))
{Undefined function or variable 'rows'.
} 
size(unique(GEOMETRY,'rows'))

ans =

   100     2

aa=unique(GEOMETRY,'rows');
format short
[GEOMETRY aaa]
{Undefined function or variable 'aaa'.
} 
[GEOMETRY aa]

ans =

    3.0000         0    1.0005   -0.0317
    2.9980    0.0634    1.0005    0.0317
    2.9920    0.1266    1.0045   -0.0951
    2.9819    0.1893    1.0045    0.0951
    2.9679    0.2511    1.0126   -0.1580
    2.9501    0.3120    1.0126    0.1580
    2.9284    0.3717    1.0246   -0.2203
    2.9029    0.4298    1.0246    0.2203
    2.8738    0.4862    1.0405   -0.2817
    2.8413    0.5406    1.0405    0.2817
    2.8053    0.5929    1.0603   -0.3420
    2.7660    0.6428    1.0603    0.3420
    2.7237    0.6901    1.0839   -0.4009
    2.6785    0.7346    1.0839    0.4009
    2.6306    0.7761    1.1112   -0.4582
    2.5801    0.8146    1.1112    0.4582
    2.5272    0.8497    1.1420   -0.5137
    2.4723    0.8815    1.1420    0.5137
    2.4154    0.9096    1.1763   -0.5671
    2.3569    0.9341    1.1763    0.5671
    2.2969    0.9549    1.2139   -0.6182
    2.2358    0.9718    1.2139    0.6182
    2.1736    0.9848    1.2547   -0.6668
    2.1108    0.9938    1.2547    0.6668
    2.0476    0.9989    1.2985   -0.7127
    1.9841    0.9999    1.2985    0.7127
    1.9208    0.9969    1.3451    0.7557
    1.8577    0.9898    1.3451   -0.7557
    1.7952    0.9788    1.3944   -0.7958
    1.7335    0.9638    1.3944    0.7958
    1.6729    0.9450    1.4461   -0.8326
    1.6137    0.9224    1.4461    0.8326
    1.5559    0.8960    1.5000   -0.8660
    1.5000    0.8660    1.5000    0.8660
    1.4461    0.8326    1.5559    0.8960
    1.3944    0.7958    1.5559   -0.8960
    1.3451    0.7557    1.6137   -0.9224
    1.2985    0.7127    1.6137    0.9224
    1.2547    0.6668    1.6729   -0.9450
    1.2139    0.6182    1.6729    0.9450
    1.1763    0.5671    1.7335   -0.9638
    1.1420    0.5137    1.7335    0.9638
    1.1112    0.4582    1.7952   -0.9788
    1.0839    0.4009    1.7952    0.9788
    1.0603    0.3420    1.8577   -0.9898
    1.0405    0.2817    1.8577    0.9898
    1.0246    0.2203    1.9208   -0.9969
    1.0126    0.1580    1.9208    0.9969
    1.0045    0.0951    1.9841    0.9999
    1.0005    0.0317    1.9841   -0.9999
    1.0005   -0.0317    2.0476   -0.9989
    1.0045   -0.0951    2.0476    0.9989
    1.0126   -0.1580    2.1108   -0.9938
    1.0246   -0.2203    2.1108    0.9938
    1.0405   -0.2817    2.1736   -0.9848
    1.0603   -0.3420    2.1736    0.9848
    1.0839   -0.4009    2.2358   -0.9718
    1.1112   -0.4582    2.2358    0.9718
    1.1420   -0.5137    2.2969   -0.9549
    1.1763   -0.5671    2.2969    0.9549
    1.2139   -0.6182    2.3569   -0.9341
    1.2547   -0.6668    2.3569    0.9341
    1.2985   -0.7127    2.4154   -0.9096
    1.3451   -0.7557    2.4154    0.9096
    1.3944   -0.7958    2.4723   -0.8815
    1.4461   -0.8326    2.4723    0.8815
    1.5000   -0.8660    2.5272   -0.8497
    1.5559   -0.8960    2.5272    0.8497
    1.6137   -0.9224    2.5801   -0.8146
    1.6729   -0.9450    2.5801    0.8146
    1.7335   -0.9638    2.6306   -0.7761
    1.7952   -0.9788    2.6306    0.7761
    1.8577   -0.9898    2.6785   -0.7346
    1.9208   -0.9969    2.6785    0.7346
    1.9841   -0.9999    2.7237   -0.6901
    2.0476   -0.9989    2.7237    0.6901
    2.1108   -0.9938    2.7660   -0.6428
    2.1736   -0.9848    2.7660    0.6428
    2.2358   -0.9718    2.8053   -0.5929
    2.2969   -0.9549    2.8053    0.5929
    2.3569   -0.9341    2.8413   -0.5406
    2.4154   -0.9096    2.8413    0.5406
    2.4723   -0.8815    2.8738   -0.4862
    2.5272   -0.8497    2.8738    0.4862
    2.5801   -0.8146    2.9029   -0.4298
    2.6306   -0.7761    2.9029    0.4298
    2.6785   -0.7346    2.9284   -0.3717
    2.7237   -0.6901    2.9284    0.3717
    2.7660   -0.6428    2.9501   -0.3120
    2.8053   -0.5929    2.9501    0.3120
    2.8413   -0.5406    2.9679   -0.2511
    2.8738   -0.4862    2.9679    0.2511
    2.9029   -0.4298    2.9819   -0.1893
    2.9284   -0.3717    2.9819    0.1893
    2.9501   -0.3120    2.9920   -0.1266
    2.9679   -0.2511    2.9920    0.1266
    2.9819   -0.1893    2.9980   -0.0634
    2.9920   -0.1266    2.9980    0.0634
    2.9980   -0.0634    3.0000   -0.0000
    3.0000   -0.0000    3.0000         0

GEOMETRY(1,:) == GEOMETRY(end,:)

ans =

  1×2 <a href="matlab:helpPopup logical" style="font-weight:bold">logical</a> array

   1   0

GEOMETRY(1,:)

ans =

     3     0

GEOMETRY(end,:)

ans =

    3.0000   -0.0000

GEOMETRY(1,2)-GEOMETRY(end,2)

ans =

   2.4493e-16

eps

ans =

   2.2204e-16

abs(GEOMETRY(1,1)-GEOMETRY(end,1)) < 1E-10 && ...
     abs(GEOMETRY(1,2)-GEOMETRY(end,2)) < 1E-10

ans =

  <a href="matlab:helpPopup logical" style="font-weight:bold">logical</a>

   1

if abs(GEOMETRY(1,1)-GEOMETRY(end,1)) < 1E-10 && ...
     abs(GEOMETRY(1,2)-GEOMETRY(end,2)) < 1E-10
    GEOMETRY=GEOMETRY(1:end-1,:);
end
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
if system_dependent('IsDebugMode')==1, dbquit; end
main_MESH
dos(['gmsh -2 ' filename '.geo']);
Info    : Running 'gmsh -2 GMSH_GEOMETRY.geo' [Gmsh 2.16.0, 1 node, max. 1 thread]
Info    : Started on Mon Feb 26 15:52:56 2018
Info    : Reading 'GMSH_GEOMETRY.geo'...
Info    : Done reading 'GMSH_GEOMETRY.geo'
Info    : Finalized high order topology of periodic connections
Info    : Meshing 1D...
Info    : Meshing curve 1 (Line)
Info    : Meshing curve 2 (Line)
Info    : Meshing curve 3 (Line)
Info    : Meshing curve 4 (Line)
Info    : Meshing curve 5 (Line)
Info    : Meshing curve 6 (Line)
Info    : Meshing curve 7 (Line)
Info    : Meshing curve 8 (Line)
Info    : Meshing curve 9 (Line)
Info    : Meshing curve 10 (Line)
Info    : Meshing curve 11 (Line)
Info    : Meshing curve 12 (Line)
Info    : Meshing curve 13 (Line)
Info    : Meshing curve 14 (Line)
Info    : Meshing curve 15 (Line)
Info    : Meshing curve 16 (Line)
Info    : Meshing curve 17 (Line)
Info    : Meshing curve 18 (Line)
Info    : Meshing curve 19 (Line)
Info    : Meshing curve 20 (Line)
Info    : Meshing curve 21 (Line)
Info    : Meshing curve 22 (Line)
Info    : Meshing curve 23 (Line)
Info    : Meshing curve 24 (Line)
Info    : Meshing curve 25 (Line)
Info    : Meshing curve 26 (Line)
Info    : Meshing curve 27 (Line)
Info    : Meshing curve 28 (Line)
Info    : Meshing curve 29 (Line)
Info    : Meshing curve 30 (Line)
Info    : Meshing curve 31 (Line)
Info    : Meshing curve 32 (Line)
Info    : Meshing curve 33 (Line)
Info    : Meshing curve 34 (Line)
Info    : Meshing curve 35 (Line)
Info    : Meshing curve 36 (Line)
Info    : Meshing curve 37 (Line)
Info    : Meshing curve 38 (Line)
Info    : Meshing curve 39 (Line)
Info    : Meshing curve 40 (Line)
Info    : Meshing curve 41 (Line)
Info    : Meshing curve 42 (Line)
Info    : Meshing curve 43 (Line)
Info    : Meshing curve 44 (Line)
Info    : Meshing curve 45 (Line)
Info    : Meshing curve 46 (Line)
Info    : Meshing curve 47 (Line)
Info    : Meshing curve 48 (Line)
Info    : Meshing curve 49 (Line)
Info    : Meshing curve 50 (Line)
Info    : Meshing curve 51 (Line)
Info    : Meshing curve 52 (Line)
Info    : Meshing curve 53 (Line)
Info    : Meshing curve 54 (Line)
Info    : Meshing curve 55 (Line)
Info    : Meshing curve 56 (Line)
Info    : Meshing curve 57 (Line)
Info    : Meshing curve 58 (Line)
Info    : Meshing curve 59 (Line)
Info    : Meshing curve 60 (Line)
Info    : Meshing curve 61 (Line)
Info    : Meshing curve 62 (Line)
Info    : Meshing curve 63 (Line)
Info    : Meshing curve 64 (Line)
Info    : Meshing curve 65 (Line)
Info    : Meshing curve 66 (Line)
Info    : Meshing curve 67 (Line)
Info    : Meshing curve 68 (Line)
Info    : Meshing curve 69 (Line)
Info    : Meshing curve 70 (Line)
Info    : Meshing curve 71 (Line)
Info    : Meshing curve 72 (Line)
Info    : Meshing curve 73 (Line)
Info    : Meshing curve 74 (Line)
Info    : Meshing curve 75 (Line)
Info    : Meshing curve 76 (Line)
Info    : Meshing curve 77 (Line)
Info    : Meshing curve 78 (Line)
Info    : Meshing curve 79 (Line)
Info    : Meshing curve 80 (Line)
Info    : Meshing curve 81 (Line)
Info    : Meshing curve 82 (Line)
Info    : Meshing curve 83 (Line)
Info    : Meshing curve 84 (Line)
Info    : Meshing curve 85 (Line)
Info    : Meshing curve 86 (Line)
Info    : Meshing curve 87 (Line)
Info    : Meshing curve 88 (Line)
Info    : Meshing curve 89 (Line)
Info    : Meshing curve 90 (Line)
Info    : Meshing curve 91 (Line)
Info    : Meshing curve 92 (Line)
Info    : Meshing curve 93 (Line)
Info    : Meshing curve 94 (Line)
Info    : Meshing curve 95 (Line)
Info    : Meshing curve 96 (Line)
Info    : Meshing curve 97 (Line)
Info    : Meshing curve 98 (Line)
Info    : Meshing curve 99 (Line)
Info    : Done meshing 1D (0.03125 s)
Info    : Meshing 2D...
Info    : Meshing surface 1 (Plane, Delaunay)
Info    : Done meshing 2D (0.0616646 s)
Info    : 1077 vertices 2251 elements
Info    : Writing 'GMSH_GEOMETRY.msh'...
Info    : Done writing 'GMSH_GEOMETRY.msh'
Info    : Stopped on Mon Feb 26 15:52:56 2018
[mesh]=readGMSH(filename);
[p,t]=readGMSH(filename);
ntri=size(t,1);
lati=zeros(1,2);
lati_tot=zeros(3*size(t,1),2);
for i=1:ntri
    triangle=t(i,:);
    lati_loc=[triangle(1) triangle(2);
        triangle(2) triangle(3);
        triangle(3) triangle(1)];
    if i==1
        lati=lati(2:end,:);
    end
    lati_tot(3*i-2:3*i,:)=lati_loc;
end
lati_tot=sort(lati_tot,2);
[e,~,~]=unique(lati_tot,'rows');
[p,e,t]=readGMSH(filename);
help pdemesh
 <strong>pdemesh</strong> Plot a PDE mesh.
    <strong>pdemesh</strong>(P,E,T) plots the mesh specified by P, E, and T.
 
    <strong>pdemesh</strong>(P,E,T,U) plots the solution column vector U using
    a mesh plot.  If U is a column vector, node data is
    assumed. If U is a row vector, element data is assumed.
    This command plots the solution substantially faster than the
    PDESURF command.
 
    H=<strong>pdemesh</strong>(P,E,T) additionally returns handles to the plotted
    graphics objects.
 
    <strong>pdemesh</strong>(MSH,...) Allows the mesh to be defined by a pde.FEMesh.
 
    <strong>pdemesh</strong>(PDEM,...) Allows the mesh to be defined by PDEM, a 
    pde.PDEModel object that contains the mesh.
 
    <strong>pdemesh</strong>(...,'Name1',Value1, 'Name2',Value2) uses name-value pairs to 
    support the following mesh plot options when plotting the mesh in
    isolation without the solution.
    
    Name        Value/{Default}   Description
    ----------------------------------------------------------------
    NodeLabels      {off} | on    - Display node labels.
    ElementLabels   {off} | on    - Display element labels.
    FaceAlpha       data          - Transparency of faces.
                                    Scalar in the range [0, 1], default = 1
 
 
    See also <a href="matlab:help initmesh">initmesh</a>, <a href="matlab:help createpde">createpde</a>, <a href="matlab:help pde.PDEModel/generateMesh">pde.PDEModel/generateMesh</a>, <a href="matlab:help pde.FEMesh">pde.FEMesh</a>
             <a href="matlab:help pdeplot">pdeplot</a>, <a href="matlab:help pdeplot3D">pdeplot3D</a>

    <a href="matlab:doc pdemesh">Reference page for pdemesh</a>

figure
pdemesh(p,e,t)
{Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('pdeplot3D', 'C:\Program Files\MATLAB\R2017a\toolbox\pde\pdeplot3D.m', 76)" style="font-weight:bold">pdeplot3D</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\pde\pdeplot3D.m',76,0)">line 76</a>)
Param value pairs expected.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('pdemesh', 'C:\Program Files\MATLAB\R2017a\toolbox\pde\pdemesh.m', 102)" style="font-weight:bold">pdemesh</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\pde\pdemesh.m',102,0)">line 102</a>)
     hh=pdeplot3D(argsToPass{:});
} 
pdemesh(p',e',t')
{Index exceeds matrix dimensions.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('pdeplot', 'C:\Program Files\MATLAB\R2017a\toolbox\pde\pdeplot.m', 269)" style="font-weight:bold">pdeplot</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\pde\pdeplot.m',269,0)">line 269</a>)
        T=sparse(t([1 4 2 5 3 6],:),t([4 2 5 3 6
        1],:),1,np,np);

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('pdemesh', 'C:\Program Files\MATLAB\R2017a\toolbox\pde\pdemesh.m', 100)" style="font-weight:bold">pdemesh</a> (<a href="matlab: opentoline('C:\Program Files\MATLAB\R2017a\toolbox\pde\pdemesh.m',100,0)">line 100</a>)
     hh=pdeplot(argsToPass{:});
} 
help triplot
 <strong>triplot</strong> Plots a 2D triangulation
    <strong>triplot</strong>(TRI,X,Y) displays the triangles defined in the
    M-by-3 matrix TRI.  A row of TRI contains indices into X,Y that
    define a single triangle. The default line color is blue.
 
    <strong>triplot</strong>(TR) displays the triangles in the triangulation TR.
 
    <strong>triplot</strong>(...,COLOR) uses the string COLOR as the line color.
 
    H = <strong>triplot</strong>(...) returns a line handle representing the displayed
    triangles edges.
 
    <strong>triplot</strong>(...,'param','value','param','value'...) allows additional
    line param/value pairs to be used when creating the plot.
 
    Example 1:
        X = rand(10,2);
        dt = delaunayTriangulation(X);
        triplot(dt)
 
    Example 2:
        % Plotting a Delaunay triangulation in face-vertex format
        X = rand(10,2);
        dt = delaunayTriangulation(X);
        tri = dt(:,:);
        triplot(tri, X(:,1), X(:,2));
 
    See also <a href="matlab:help trisurf">trisurf</a>, <a href="matlab:help trimesh">trimesh</a>, <a href="matlab:help delaunay">delaunay</a>, <a href="matlab:help triangulation">triangulation</a>, <a href="matlab:help delaunayTriangulation">delaunayTriangulation</a>.

    <a href="matlab:doc triplot">Reference page for triplot</a>

figure
triplot(t,p(:,1),p(:,2));
figure; hold on, axis equal;
triplot(t,p(:,1),p(:,2));
