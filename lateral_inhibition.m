%{
Will Stonehouse
Psych 186B: Neural Networks
Homework 4: Lateral Inhibition
Assignment: Replicate results in Chapter 4 of Anderson text
%}

clear;
close all;
clc;

%% Figure 4.19

% Initialize parameters
Dimension = 80;
Number_of_iterations = 1000;

Lower_limit = 0;
Upper_limit = 60;

Length_const = 2;
Max_strength = 0.1;
Epsilon = 0.01;

% Computing lateral inhibition
[old, new] = lateral_inhibition(Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon);

% Plot figure
plot(new(12:30), 'k*');
hold on;
plot(old(12:30), 'k+');
hold off;

%% Figure 4.20

% Initialize parameters
Dimension = 80;
Number_of_iterations = 100;

Lower_limit = 0;
Upper_limit = 60;

Length_const = 2;
Max_strength = 0.2;
Epsilon = 0.01;

% Computing lateral inhibition
[old, new] = lateral_inhibition(Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon);

% Plot figure
plot(new(12:30), 'k*');
hold on;
plot(old(12:30), 'k+');
hold off;

%% Figure 4.21

% Initialize parameters
Dimension = 80;
Number_of_iterations = 100;

Lower_limit = 0;
Upper_limit = 60;

Length_const = 2;
Max_strength = 0.5;
Epsilon = 0.01;

% Computing lateral inhibition
[old, new] = lateral_inhibition(Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon);

% Plot figure
plot(new(12:30), 'k*');
hold on;
plot(old(12:30), 'k+');
hold off;

%% Figure 4.22

% Initialize parameters
Dimension = 80;
Number_of_iterations = 100;

Lower_limit = 0;
Upper_limit = 60;

Length_const = 2;
Max_strength = 1;
Epsilon = 0.01;

% Computing lateral inhibition
[old, new] = lateral_inhibition(Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon);

% Plot figure
plot(new(12:30), 'k*');
hold on;
plot(old(12:30), 'k+');
hold off;

%% Figure 4.23

% Initialize parameters
Dimension = 80;
Number_of_iterations = 100;

Lower_limit = 0;
Upper_limit = 60;

Length_const = 2;
Max_strength = 2;
Epsilon = 0.01;

% Computing lateral inhibition
[old, new] = lateral_inhibition(Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon);

% Plot figure
plot(new(12:30), 'k*');
hold on;
plot(old(12:30), 'k+');
hold off;

%% Figure 4.26

% Initialize parameters
Dimension = 80;
Number_of_iterations = 100;

Lower_limit = 0;
Upper_limit = 60;

Length_const = 10;
Max_strength = 1;
Epsilon = 0.05;

% Initialize state vector
Initial_state_vector = zeros(Dimension,1);
Initial_state_vector(16) = 10;
Initial_state_vector(17) = 20;
Initial_state_vector(18) = 30;
Initial_state_vector(19) = 40;
Initial_state_vector(20) = 50;
Initial_state_vector(21) = 40;
Initial_state_vector(22) = 30;
Initial_state_vector(23) = 20;
Initial_state_vector(24) = 10;

% Computing lateral inhibition
[old, new] = WTAN(Initial_state_vector, Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon);

% Plot figure
plot(new(12:30), 'k*');
hold on;
plot(old(12:30), 'k+');
hold off;

%% Figure 4.27

% Initialize parameters
Dimension = 80;
Number_of_iterations = 100;

Lower_limit = 0;
Upper_limit = 60;

Length_const = 10;
Max_strength = 1;
Epsilon = 0.055;

% Initialize state vector
Initial_state_vector = zeros(Dimension,1);
Initial_state_vector(1:16) = 10;
Initial_state_vector(17) = 20;
Initial_state_vector(18) = 30;
Initial_state_vector(19) = 40;
Initial_state_vector(20) = 50;
Initial_state_vector(21) = 40;
Initial_state_vector(22) = 30;
Initial_state_vector(23) = 20;
Initial_state_vector(24:80) = 10;

% Computing lateral inhibition
[old, new] = WTAN(Initial_state_vector, Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon);

% Plot figure
plot(new(12:30), 'k*');
hold on;
plot(old(12:30), 'k+');
hold off;

%% Figure 4.28

% Initialize parameters
Dimension = 80;
Number_of_iterations = 1000;

Lower_limit = 0;
Upper_limit = 60;

Length_const = 10;
Max_strength = 1;
Epsilon = 0.09;

% Initialize state vector
Initial_state_vector = zeros(Dimension,1);
Initial_state_vector(1:80) = 10;
Initial_state_vector(14) = 20;
Initial_state_vector(15) = 30;
Initial_state_vector(16) = 20;
Initial_state_vector(19) = 20;
Initial_state_vector(20) = 30;
Initial_state_vector(21) = 40;
Initial_state_vector(22) = 30;
Initial_state_vector(23) = 20;

% Computing lateral inhibition
[old, new] = WTAN(Initial_state_vector, Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon);

% Plot figure
plot(new(12:30), 'k*');
hold on;
plot(old(12:30), 'k+');
hold off;

%% Figure 4.29

% Initialize parameters
Dimension = 80;
Number_of_iterations = 1000;

Lower_limit = 0;
Upper_limit = 60;

Length_const = 10;
Max_strength = 2;
Epsilon = 0.1423;

% Initialize state vector
Initial_state_vector = zeros(Dimension,1);
Initial_state_vector(1:80) = 10;
Initial_state_vector(14) = 20;
Initial_state_vector(15) = 30;
Initial_state_vector(16) = 20;
Initial_state_vector(19) = 20;
Initial_state_vector(20) = 30;
Initial_state_vector(21) = 40;
Initial_state_vector(22) = 30;
Initial_state_vector(23) = 20;

% Computing lateral inhibition
[old, new] = WTAN(Initial_state_vector, Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon);

% Plot figure
plot(new(12:30), 'k*');
hold on;
plot(old(12:30), 'k+');
hold off;

%% Lateral inhibition function

function [Initial_state_vector, New_vector] = lateral_inhibition(Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon)

    % Initialize state vector
    Initial_state_vector = zeros(Dimension,1);
    Initial_state_vector(1:20) = 10;
    Initial_state_vector(21:60) = 40;
    Initial_state_vector(61:80) = 10;
    State_vector = Initial_state_vector;
    
    % Setting inhibitory weights
    Inhibitory_weights = zeros(Dimension,Dimension);
    for i = 1:Dimension
        for j = 1:Dimension
            Distance = abs(i-j);  
            if Distance <= Dimension/2
                Inhibitory_weights(i,j) = - Max_strength * exp(-Distance/Length_const);
            end
            if Distance > Dimension/2
                Distance = Dimension-Distance;
                Inhibitory_weights(i,j) = - Max_strength * exp(-Distance/Length_const);
            end
        end
    end
    
    % Computing final state vector
    for i=1:Number_of_iterations
        for j=1:Dimension
            new = 0;
            for k=1:Dimension
                new = new + Inhibitory_weights(j,k) * State_vector(k);
            end
            result = State_vector(j) + Epsilon * (Initial_state_vector(j) + new - State_vector(j));
            if result < Lower_limit
                result = Lower_limit;

            elseif result > Upper_limit
                result = Upper_limit;
            end
            New_vector(j) = result;
        end
        State_vector = New_vector;
    end
end

%% Winner Takes All Networks

function [Initial_state_vector, New_vector] = WTAN(Initial_state_vector, Dimension, Number_of_iterations, Lower_limit, Upper_limit, Length_const, Max_strength, Epsilon)
    
    State_vector = Initial_state_vector;
    
    % Setting inhibitory weights
    Inhibitory_weights = zeros(Dimension,Dimension);
    for i = 1:Dimension
        for j = 1:Dimension
            Distance = abs(i-j);  
            if Distance <= Dimension/2
                if i==j
                    Inhibitory_weights(i,j) = 0;
                else
                    Inhibitory_weights(i,j) = - Max_strength * exp(-Distance/Length_const);
                end
            end
            if Distance > Dimension/2
                if i==j
                    Inhibitory_weights(i,j) = 0;
                else
                    Distance = Dimension-Distance;
                    Inhibitory_weights(i,j) = - Max_strength * exp(-Distance/Length_const);
                end
            end
        end
    end
    
    % Computing final state vector
    for i=1:Number_of_iterations
        for j=1:Dimension
            new = 0;
            for k=1:Dimension
                new = new + Inhibitory_weights(j,k) * State_vector(k);
            end
            result = State_vector(j) + Epsilon * (Initial_state_vector(j) + new - State_vector(j));
            if result < Lower_limit
                result = Lower_limit;
            elseif result > Upper_limit
                result = Upper_limit;
            end
            New_vector(j) = result;
        end
        State_vector = New_vector;
    end
end