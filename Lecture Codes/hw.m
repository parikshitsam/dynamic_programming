%while loop hw example
clc
clear
bal = 5000;
year = 0;
while bal<1000000;
    bal=1.08*bal + 5000;
    year= year+1;
end
 disp(year)

%if loop hw example
clc
clear
max_coffee=input('how many cups of coffee you had today?');


     if(max_coffee>=3)
         fprintf('you had max cups today')
     elseif(max_coffee<3 & max_coffee>=2)
         fprintf('you have one cup left')
     else
         fprintf('you have two cups left')
     end
     
     