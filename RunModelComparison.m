function [sumprobability1,sumprobability2]=RunModelComparison(iterations)
% This function computes the probability that mutational changes produce
% the desired phenotypic effect in two different pathways. The input
% variable "iterations" is the number of samples of the reaction kinetics.
% More "iterations" will produce better estimates of the probabilities.
% The output of the function is two matrices, one for each pathway, that
% contain the probability estimate for different combinations of the
% probability of enabling (rows) and disabling (columns) change.

perange=10.^([-6:.1:-2]');
pdrange=perange;
sumprobability1=zeros(length(perange),length(pdrange));
for i1=1:iterations
    sampleprobability=PhenotypeChange1(perange,pdrange);
    sumprobability1=sumprobability1+sampleprobability;
end

sumprobability2=zeros(length(perange),length(pdrange));
for i1=1:iterations
    sampleprobability=PhenotypeChange2(perange,pdrange);
    sumprobability2=sumprobability2+sampleprobability;
end
end


function totalprob=PhenotypeChange1(perange,pdrange)
    % This is a helper function of RunModelComparison that computes the
    % probability that mutational changes produce the desired phenotypic
    % effect for a range of probabilities of enabling and disabling change.
    % For illustration we use the AWS pathway's production of wrinkly
    % spreader phenotypes
    concentrations=SampleAWSConcentrations();
    rates=SampleAWSRates();
    time=SampleTime();
    indicator=SolveAWSModel(concentrations,rates,time);
    numreactions=length(rates);
    totalprob=0;
    for k=1:3^numreactions
        mutations=NumberToMutationArray(numreactions,k);
        newrates=Mutate(rates,mutations);
        newindicator=SolveAWSModel(concentrations,newrates,time);
        if Changed(indicator,newindicator)
            eventprob=ComputeProbability(mutations,perange,pdrange);
            totalprob=totalprob+eventprob;
        end
    end
end

function eventprob=ComputeProbability(mutations,probenhancerange,probdeleterange)
    % This is a helper function for PhenotypeChange1 and PhenotypeChange2 that calculates the probability that a specific set of mutations will occur.
    eventprob=ones(length(probenhancerange),length(probdeleterange));
    numenhance=sum(mutations==1);
    eventprob=eventprob.*(probenhancerange.^numenhance*ones(1,length(probdeleterange)));    
    numdelete=sum(mutations==-1);
    eventprob=eventprob.*(ones(length(probenhancerange),1)*(probdeleterange').^numdelete);    
    numnothing=sum(mutations==0);
    probnothing=ones(length(probenhancerange),length(probdeleterange))-ones(length(probenhancerange),1)*(probdeleterange')-probenhancerange*ones(1,length(probdeleterange));
    eventprob=(probnothing.^numnothing).*eventprob;
end

function AWSrates=SampleAWSRates()
    % This is a helper function for PhenotypicChange1 that samples the
    % reaction rates for pathway1.
    AWSrates=10.^(4*(rand(4,1)-.5));
end


function time=SampleTime()
    % This is a helper function for PhenotypicChange1 that samples the time for the differential equation solver.
    time=10;
end


function AWSconcentrations=SampleAWSConcentrations()
    % This is a helper function for PhenotypicChange1 that samples the
    % concentrations for pathway1.
    AWSconcentrations=rand(8,1)*10;
end




function mutations=NumberToMutationArray(numreactions,num)
    % This is a helper function for PhenotypeChange1 and PhenotypeChange2 that converts an integer "num" into a specific set of enabling/disabling/no change mutations in "numreactions" reactions.
	mutations=zeros(numreactions,1);
	for i1=numreactions-1:-1:0;
		mutations(i1+1)=floor(num/(3^i1));
		num=num-mutations(i1+1)*3^i1;
	end
	mutations=mutations-1;
end


function AWSindicator=SolveAWSModel(AWSconcentrations,AWSrates,time)
    % This is a helper function for PhenotypicChange1 that solves the
    % dynamical system and returns the indicator concentration
    [timevalues,solution]=ode45(@(t,y) odeAWS(t,y,AWSrates),[0 time],AWSconcentrations);
    AWSindicator=solution(end,7);
end

function yp=odeAWS(t,y,AWSrates)
    % This is a helper function for AWSindicator that describes the
    % dynamical system of pathway 1 that will be solved numerically. For
    % illustration we write the AWS pathway for generating a wrinkly
    % spreader phenotype.
    X=y(1);
    XR=y(2);
    O=y(3);
    Op=y(4);
    OX=y(5);
    R=y(6);
    RR=y(7);
    Si=y(8);
    r1=AWSrates(1);r2=AWSrates(2);
    r3=AWSrates(3);r4=AWSrates(4);
    yp=zeros(size(y));
    yp(1)=-r2*X*Op-r3*X*R; % dX/dt
    yp(2)=r3*X*R; % dXR/dt
    yp(3)=-r1*Si*O; % dO/dt
    yp(4)=r1*Si*O-r2*X*Op; % dOp/dt
    yp(5)=r2*Op*X; % dOX/dt
    yp(6)=-r3*X*R-r4*R*R; % dR/dt
    yp(7)=r4*R*R; % dRR/dt
    yp(8)=0;
end

function verdict = Changed(indicator,newindicator)
    % This is a helper function for PhenotypeChange1 and possibly
    % PhenotypeChange2 that returns a boolean variable signifying if the
    % change in indicator's concentration is enough to trigger a phenotypic
    % change.
    if newindicator-indicator>.000001
        verdict=true;
    else
        verdict=false;
    end
end

function newrates=Mutate(rates,mutations)
    % This is a helper function for PhenotypeChange1 and PhenotypeChange2
    % that takes in a vector of reaction rates "rates" and a set of
    % mutations "mutations" and produces a vector for the new, altered
    % reaction rates "newrates".
    factors=[10.^(-2*rand()) 1 10.^(2*rand())];
    newrates=zeros(size(rates));
    for i1=1:length(rates)
        switch mutations(i1)
            case -1
                newrates(i1)=10.^(-2*rand())*rates(i1);
            case 0
                newrates(i1)=rates(i1);
            case 1
                newrates(i1)=10.^(2*rand())*rates(i1);
        end
    end
end

function totalprob=PhenotypeChange2(perange,pdrange)
    % This function should mirror PhenotypeChange1 except describing
    % another pathway. Here, we simply use PhenotypeChange1 in order to
    % generate output for this sample code. Note the contour plot generated
    % by PlotModelComparison should thus show a similar likelihood for the two
    % pathways.
    totalprob=PhenotypeChange1(perange,pdrange);
end

