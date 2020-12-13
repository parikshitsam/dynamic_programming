%input command
clc
clear
y=input('enter first value of y = '); %it will ask for input value
p=input('enter first value of p = '); %it will ask for input value
z=y+p;
disp(z) %display the result in command window


%fprintf
clc
clear
fprintf('i am studying macro') %print/display result on command window
fprintf('today is friday,i am enjoying attending ds')
fprintf('10292863')


%while loop example
clc
clear
N=0;  %intializing value 
sum_x=0; %intializing value
x=input('enter first value;'); %it will ask for input value

while x>=0  %specifying condition
    sum_x = sum_x + x;
    N = N+1;
    x = input('enter next value;');
 
end
mean_x = sum_x/N; %calculating mean 
fprintf('the mean of entered data is: %f', mean_x); %displaying mean in command window 


%for loop example-1
clc
clear
 for i=1:10 %I want to run this loop 10 times
     fprintf('i am enjoying coding today\n'); %I want each result to be printed in a new line
 end
 
%for loop example-2 (calculating factorial)
clc
clear 
n= 70; %value of n
fact = 1; %intializing value

for i =1:n
    fact = fact*i;
end
fprintf('factorial is: %f', fact);


%if condition-example 1
clc
clear
marks=input('enter the percentage marks:'); %asking students to enter marks
if (marks>=30) %specifying condition
    fprintf('pass')
else
    fprintf('fail');
end
 

%if condition-example 2
clc
clear
marks=input('enter the percentage marks:');

if(marks>=60)
    fprintf('first division')
elseif (marks>=45 && marks<60)
    fprintf('second division')
elseif (marks>=30 && marks<45)
    fprintf('third division')
else
    fprintf('fail');
end


%if under if condition
clc
clear
marks=input('enter the percentage marks:');
if(marks>=60)
    fprintf('first division')
else
     if(marks>=45 && marks<60)
    fprintf('second division')
     else
    if(marks>=30 && marks<45)
    fprintf('third division')
    else
    if (marks<30)
        fprintf ('fail');
    end
    end
    end
end

