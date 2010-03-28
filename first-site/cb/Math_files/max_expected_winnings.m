
function max_expected_winnings(moneystart,increment,p,rounds,nprob,decreasement,moneyorgame)

% This function calculates the maximum expected winnings or 
% game in whic the maximum winnings were earned for a standard
% gamblers ruin problem where the probability of winning is 
% varied from its initial amount p, to observe it sensitivity
%
% author: Cameron Bracken, 04-05-07
%
% input variable list:
% moneystart=         amout of money you start with 
% increment=          amount of money you win or lose each game (fixed)
% p=                  initial probability of sucess, must be <.5
% rounds=             number of rounds to average over
% nprob=              number of probilities to plot with corresponding game
%                     or max winnings
% decreasement=       amount to decrease p after one simulation 
% moneyorgame=        character, 'money' will plot the max expected
%                     winnings and 'game' will plot expected game in 
%                     which the max money was
% earned
% 
% This is a sample input:
% max_expected_winnings(100,1,.49,1000,40,.01,'money')

hold on

out=zeros(nprob,3);

for j=1:nprob
    maxmoneyave=0;
    maxgameave=0;
    probcount=j;
    for i=1:rounds 
        money=moneystart;
        maxmoney=money;
        game=0;
        maxgame=1;
        while money>0
            r=rand(1);
            if r<=p
                money=money+increment;
                game=game+1;
                if money>maxmoney
                    maxmoney=money;
                    maxgame=game;
                end
            else
                money=money-increment;
                game=game+1;
            end
    %         if money<=0
    %             y=['busted on game...max winnings were...max winnings on game'] 
    %             x=[game max maxgame]        
    %             break
    %         end 
        end
        maxmoneyave=maxmoneyave+maxmoney;
        maxgameave=maxgameave+maxgame;
    end

    maxmoneyave=maxmoneyave/rounds;
    maxgameave=maxgameave/rounds;
    
    out(j,:)=[p maxmoneyave maxgameave];
    p=p-decreasement;
    if p>0
        p
    else 
        break
    end
end

if strcmp(moneyorgame,'money')
    plot(out(1:probcount,1),out(1:probcount,2),'k')
    %axis([0 0.5 moneystart max(out(:,2))]);
    xlabel('p')
    ylabel(['Average winnings ($)'])
    text(.1,moneystart+1,['Starting Money $',num2str(moneystart)])
    text(.1,moneystart+2,['Increment $',num2str(increment)])
elseif strcmp(moneyorgame,'game')
    plot(out(1:probcount,1),out(1:probcount,3),'k')
    %axis([0 0.5 moneystart max(out(:,2))]);
    xlabel('p')
    ylabel(['Average # of games to max money'])
    text(.1,10,['Starting Money $',num2str(moneystart)])
    text(.1,20,['Increment $',num2str(increment)])
end

