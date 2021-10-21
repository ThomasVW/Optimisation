%%%%%%%%%%%%% Projet d’optimisation continue %%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Estimation robuste: les M-estimateurs %%%%%%%%%%%%%
%% 1) Description du problème
clear;
close all; 

file = open('data.mat');

x = file.x;
y= file.y_noisy;

%% 2) Estimation au sens des moindres carrés

plot(x,y,'ro');
title('Data points');

%Q1 On chosis défénit une plage de valeurs
taille=1000;
a = linspace(-10,20,taille);
b = linspace(-20,20,taille);
C_moindres_carres = @(a_, b_) sum( (a_ * x + b_ - y) .* (a_ * x + b_ - y) );

C = ones(length(x));
for i=1:1:taille
    for j=1:1:taille
        C(i,j) = C_moindres_carres(a(j), b(i));
    end
end

subplot(2,2,1);
mesh(a,b,C);
title('mesh');
colorbar;
subplot(2,2,2);
contour(a,b,C,50);
title('contour');
%%
% Q2 
% estimation de l'argmin de a et b  -> on part d'un echantillonage entre
% -50 et 50 pour (a,b) et on diminue l'intervalle jusqu' à trouver une
% estimation des paramètres a et b, par echantillonnage régulier on obtient
% a€[5,8] et b€[-4,-1]
%%
%Q3
%somme de 1 à n des (a*xk + b - yk)^2

X = ones(size(x,1), 2);
X(:,1) = x;
Y = y;

AB = ( inv(X.' * X) ) * X.' * Y;

hold on,
plot(AB(1), AB(2), '+r'),
hold off;

%%
%Q4 et Q5
[M,I] = min(C(:));
[I_row, I_col] = ind2sub(size(C),I);
a(I_row)
b(I_col)

% affichage de la courbe résultante

y_estim = a(I_row)*x + b(I_col);
subplot(2,2,3);
plot(x,y_estim,'y--', x,y,'ro');
title('Optimisation');
legend('ajustement au sens des moindres carrés','points de mesure');


%% Estimation robuste
%On calcule la matrice gradient à partir des données de la formule de l énnoncé.

%% Question 6
% Nous allons utiliser la méthode des plus fortes pente avec
% préconditionnement. 
% Pour cela, le pas va être déterminé de manière linéaire à chaque
% itération de façon à respecter les conditions de Wolf.
% Le pas ne doit pas être trop grand, ni trop petit. Pour cela, nous allons
% utilisé l'algortihme de Fletcher et Lemaéchal.

% dataset : x,y
% Fonction de coût : C(a,b)
% Gradient : G(a,b,:)

%% Q06

ab_k = [-10; 10];

epsilon = 1e-6;
list_ab = zeros(1000, 2);

list_ab(1, :) = ab_k;
k = 2;

f = @(x_k) norm( (X * x_k - Y) )^2 ;
grad_f = @(x_k) 2 * X.' *( X * x_k - Y );

while ( norm( grad_f(ab_k) ) > epsilon )
    d_k = - grad_f(ab_k);
    alpha_k = Fletcher_Lemarechal(f, grad_f, ab_k); %Algo fletcher
    ab_k = ab_k + alpha_k * d_k;
    list_ab(k, :) = ab_k;
    k = k + 1;
end

figure(2),
hold on,
contour (a, b , C, 50), 
colorbar, 
xlabel('a'), 
ylabel('b'),
title('Recherche du minimum'),
plot(AB(1), AB(2), '+r'),
plot(list_ab(1:k-1, 1), list_ab(1:k-1, 2)),
legend('C moindres carres', 'Q2', 'Fletcher Lemarechal'),
hold off;

%% Q07

ro = @(r, sigma) 0.5 * log( 1 + (r./sigma).^2 );
d_ro = @(r, sigma) r./(sigma^2+r.^2);
d2_ro = @(r, sigma) (sigma^2-r.^2)/(sigma^2+(r.^2)).^2;
r = linspace(-10, 10, 1000);

figure(4),
hold on,
plot(r, ro(r,1),'b'),
plot(r, d_ro(r,1),'r'),
plot(r, d2_ro(r,1),'y'),
axis([-10 10 -0.5 3]);
hold off;


%% Q08
robuste = zeros(taille);

for i=1:taille
  for j=1:taille
    robuste(i, j) = PenalisationDeCauchy(1, a(i), b(j), x, y);
  end
end
[M,I] = min(robuste(:));
[I_row, I_col] = ind2sub(size(robuste),I);
a(I_row)
b(I_col)

figure(7);
contour(a, b, robuste, 40); xlabel('b'); ylabel('a');
legend({'Fonction coût'});
title('Modélisation de l''estimation robuste');

% Q09
aGradient = zeros(taille);
bGradient = zeros(taille);

for i=1:taille
  for j=1:taille
    aGradient(i, j) = AGradient(1, a(i), b(j), x, y);
    bGradient(i, j) = BGradient(1, a(i), b(j), x, y);
  end
end
%
 %On crée deux vecteurs qui vont permettrent de servir de repère à la fonction 
 %quiver. Ils vont balayer tous les points de l'espace en y référencant le 
 %gradient correspondant.
redimensionnerMatrix = zeros(1, taille^2);
redimensionnerMatrix2 = zeros(1, taille^2);

for i=1:taille
  redimensionnerMatrix(taille*(i-1)+1:taille*i) = linspace(a(i), a(i), taille);
  redimensionnerMatrix2(taille*(i-1)+1:taille*i) = a;
end

 %On met les matrices gradient à plat sous forme de vecteur.
AGrad = reshape(aGradient, [1, taille^2]);
BGrad = reshape(bGradient, [1, taille^2]);

hold on;
quiver(redimensionnerMatrix, redimensionnerMatrix2, AGrad, BGrad);
axis equal;
title('Gradient (a b) de la fonction coût');
hold off;


%% Q10

%% Q11

%% Q12

%% Q13

%% Q14

%% Q15
