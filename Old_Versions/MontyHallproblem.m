% Monty Hall problem
clear; clc

% Switch Choice
x=[1 2 3];
win=0;
for i=1:10000
x=randperm(length(x));
xnew=x;

choice=randi([1 3]); 

if choice==1
    win=win+1;
    continue
end

xnew(xnew==choice)=[];

xnew(xnew==randi(sort(xnew)))=[];

% ndchoice=randsample(xnew,1);

if xnew==1
    win=win+1;
end
end

disp(['Switch Choice and win ',num2str((win/10000)*100),' Percent of the time!'])


% Stay with choice
x=[1 2 3];
win=0;
for i=1:10000
x=randperm(length(x));
xnew=x;

choice=randi([1 3]); 

if choice==1
    win=win+1;
    continue
end
end

disp(['Keep choice and win ',num2str((win/10000)*100),' Percent of the time!'])
