%La presente simulación es para un sistema MISO que se compone de 3 LEDs
%transmisores y 1 receptor óptico, se ambienta en un túnel minero de 4
%metros de altura, 4 metros de ancho y 8 metros de largo.

%Consideraciones:
%El túnel tiene lumianrias colineales entre sí a lo largo del eje X
%(8[m]), se estudia la separación que deben tener las luminarias entre sí
%por medio de un estudio luminotécnico tomando consideraciones básicas del
%reglamento de seguridad minera y el manual de carreteras, y así también se
%define la potencia de transmisión en cada caso para no generar
%deslumbramiento.
%Los supuestos son ideales y basados en la documentación de fotodetectores
%utilizados en investigaciones similares dentro del estudio de VLC en la
%universidad Diego Portales. Sin embargo, son parámetros que pueden ser
%modificados para estudiar otros fenómenos.
%No se estudia el comportamiento del canal NLOS, ni ruido en el ambiente.


clc;
clear all;
close all;

% DEFINICIÓN DE PARÁMETROS INICIALES
pw = 0.6;                                    % coeficiente de reflexión del área reflectiva
ang_rad = 60;                                % semi-ángulo de mitad de potencia del LED [60-90]
m = 1;                                       % número Lambertiano
Ap = 0.0001;                                 % área física del receptor (1 cm^2) (igual a Adet)
eta = 1.5;                                   % índice de refracción del PD
fov = 70;                                    % field of view

% parámetros de los ángulos de inclinación y rotación
beta_i = 45;                                 % ángulo de inclinación del LED con respecto al eje z
alpha_i = 45;                                % ángulo de rotación del LED con respecto al eje x
beta_j = 0;                                  % ángulo de inclinación inicial del PD
alpha_j = 45;                                % ángulo de rotación del PD con respecto al eje x

% posición Ti
x_i = 3;
y_i = 0.5;
z_i = 3.6;

% posición del receptor
x_j = 3;
y_j = 1; 
z_j = 2.6;

% Calcular normal del LED
N_LED = [sin(beta_i * pi / 180) * cos(alpha_i * pi / 180), ...
          sin(beta_i * pi / 180) * sin(alpha_i * pi / 180), ...
          cos(beta_i * pi / 180)];

% Vector desde el LED al receptor
vector_to_receiver = [x_j - x_i, y_j - y_i, z_j - z_i];

% Normalizar el vector
norm_vector = vector_to_receiver / norm(vector_to_receiver);


% DEFINICIÓN DE PARÁMETROS DEL CÓDIGO SECUNDARIO
theta = 70;                             % Semi-ángulo de media potencia
Adet = Ap;                               % Área del detector
Ts = 1;                                  % Ganancia de un filtro óptico
index = 1.5;                             % Índice de refracción de una lente
FOV = 60 * pi / 180;                     % Campo de visión (FOV) del receptor
G_Con = (index^2) / sin(FOV);            % Ganancia de un concentrador óptico
c = 3e8;                                 % Velocidad de la luz en m/s


%posiciones de los LEDs a lo largo del eje X
LED_pos = [
    2, 4, 3.6;   % LED 1
    4, 4, 3.6;   % LED 2
    6, 4, 3.6    % LED 3
];

Pos_FD = zeros(8, 3);
for i = 1:8                         %en total son 8 puntos y solo avanzará en el eje X
    Pos_FD(i, :) = [i, 2, 2.6];     % Posiciones del fotodetector (x, y, z)
end

FD = [x_j, y_j, z_j];                       % Usar la posición del receptor
P_transmision = 13.1;                       % Potencia de transmisión de los LEDs
P_recepcion = zeros(1, size(LED_pos, 1));   % potencia recibida

% Calcular la potencia recibida desde cada LED en función de Hlos
for i = 1:size(LED_pos, 1)
    Xi = LED_pos(i, 1);
    Yi = LED_pos(i, 2);
    hi = LED_pos(i, 3);
    
    Xp = FD(1);
    Yp = FD(2);
    hp = FD(3);
    
    % Calcular la distancia y la potencia recibida
    D = sqrt((hi - hp)^2 + (Xi - Xp)^2 + (Yi - Yp)^2);
    H_LOS_dB =  Ts * G_Con * Adet * (m + 1) * ( abs((hi - hp)) )^(m + 1) / (2 * pi * D^2);
    P_recepcion(i) = (P_transmision * H_LOS_dB);
end

disp('Potencias recibidas en el fotodetector desde cada LED:');
disp(P_recepcion);

% Encontrar la potencia máxima y su índice
[P_recep_max, max_index] = max(P_recepcion);

% Encontrar todos los LEDs con la potencia máxima
max_indices = find(P_recepcion == P_recep_max);

% Extraer las coordenadas de los LEDs correspondientes
max_LED_coords = LED_pos(max_indices, :);

fprintf('La distancia correspondiente para dicha potencia es: %.2f metros\n', ...
    sqrt((max_LED_coords(1, 3) - FD(3))^2 + (max_LED_coords(1, 1) - FD(1))^2 + (max_LED_coords(1, 2) - FD(2))^2));
fprintf('Las coordenadas del LED con la potencia más alta son: (%.2f, %.2f, %.2f)\n', ...
    max_LED_coords(1, 1), max_LED_coords(1, 2), max_LED_coords(1, 3));

% Si hay más de un LED con la potencia máxima, mostrar las coordenadas adicionales
if length(max_indices) > 1
    for j = 2:length(max_indices)
        fprintf('Otro LED con la misma potencia máxima tiene coordenadas: (%.2f, %.2f, %.2f)\n', ...
            max_LED_coords(j, 1), max_LED_coords(j, 2), max_LED_coords(j, 3));
    end
    fprintf('La ubicación del dispositivo se encuentra entre los LEDs con coordenadas (%.2f, %.2f, %.2f) y (%.2f, %.2f, %.2f)\n', ...
        max_LED_coords(1, 1), max_LED_coords(1, 2), max_LED_coords(1, 3), ...
        max_LED_coords(2, 1), max_LED_coords(2, 2), max_LED_coords(2, 3));
end

fprintf('La potencia máxima recibida es: %.4f W\n', P_recep_max);

% Graficar las posiciones de los LEDs y el fotodetector
figure;
hold on;
grid on;

scatter3(LED_pos(:, 1), LED_pos(:, 2), LED_pos(:, 3), 100, 'filled');
scatter3(FD(1), FD(2), FD(3), 100, 'm', 'filled');


xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Posiciones de los LEDs y el Fotodetector');
legend('LEDs', 'Fotodetector');
view(3);





                
       % ----- Calcular ángulos del receptor para mejorar la recepción de señal -----%

beta_j = acosd(norm_vector(3));                     % Ángulo de inclinación
alpha_j = atan2d(norm_vector(2), norm_vector(1));   % Ángulo de rotación

disp(['Nuevo ángulo de inclinación del receptor: ', num2str(beta_j), '°']);
disp(['Nuevo ángulo de rotación del receptor: ', num2str(alpha_j), '°']);





disp('----------------------------------------------------------------------------------');  %
disp('------Resultados de potencia recibida en razón de tiempos en un sistema MISO------');
disp('----------------------------------------------------------------------------------');  % 







                         %----- Hlos con respuesta al impulso ------%



matriz_dist_ecuclidiana = zeros(3, 8);      %matriz de distancias ecuclidiana entre receptor en cada distancia en el eje X y cada uno de los transmisores de posición fija
Precepcion_t = zeros(8, 6);                 %en cada una de las posiciones almacena la altura del PD(h_i), la distancia del LED de distancia euclidiana más cercana al PD (d_cercana_i), los impulsos del canal de cana uno de los 3 LEDs (H_LOS_W_t_all(i, j)), t_i que corresponde a la respuesta al impluso de la distancia más cercana
P_recep = zeros(8, 3);                      % Inicializar P_recep para almacenar los resultados entre la Ptransmitida y H_LOS_W_t_all 

H_LOS_W_t_all = zeros(8, 3); % 8 posiciones para el FD y 3 LEDs

%parámetros de ruido en el receptor
k = 1.38e-23;   % Constante de Boltzmann (J/K)
T = 300;        %Temperatura en Kelvin
B = 100e6;      % Ancho de banda en Hz (10 MHz)
R = 50;         % Resistencia en ohmios
q = 1.602e-19;  % Carga del electrón en coulombs
I_b = 1e-3;     % Corriente de fotodetección en amperios (1 mA)


% Calcular ruido térmico y de disparo
thermal_noise = ((sqrt( 4 * k * T *R * B))^2) / R;  % Ruido térmico en W (valor en potencia V^2/R)
I_n = ((sqrt(2 * q * I_b * B))^2) * R ;             % Corriente de ruido de disparo (valor en potencia I^2*R)


% Inicializaar vectores de SNR
SNR_Paralelos_por_LED = zeros(8,3);
SNR_MAX_Paralelos = zeros(8,2); %se tienen 2 valores en cada vector por SNR y posición, esta info sirve para gráficas
SNR_MIN_Paralelos = zeros(8,2); %se tienen 2 valores en cada vector por SNR y posición, esta info sirve para gráficas


% Inicialización del vectores para almacenar la potencia con ruido.
Pr_ruido_termico_vec = zeros(8, 3); 
Pr_ruido_disparo_vec = zeros(8, 3); 
Pr_ruido_disparo_y_termico_vec = zeros(8, 3); 



for i = 1:8
    % 1. Altura del FD
    h_i = Pos_FD(i, 3);
    
    for j = 1:size(LED_pos, 1)
        % Distancia euclidiana desde el FD al LED
        distancia = norm((Pos_FD(i, :) - LED_pos(j, :)));
        matriz_dist_ecuclidiana(j, i) = distancia;
    end
    
    % Obtener menor distancia de fotodetector al los LEDs
    d_cercana_i = min(matriz_dist_ecuclidiana(:, i));
    
    % Calcular tiempo t_i
    t_i = d_cercana_i / c;
    
    % Calcular H_los para cada LED
    H_LOS_W_t = zeros(1, 3); % Inicializar en W

    for j = 1:size(LED_pos, 1)
        hi = Pos_FD(i, 3);                  % Altura del FD
        hp = LED_pos(j, 3);                 % Altura del LED
        D = matriz_dist_ecuclidiana(j, i);  % Distancia correspondiente
        
        % Cálculo de H_LOS en W
        H_LOS_W_t(j) = Ts * G_Con * Adet * (m + 1) * ( abs((hi - hp)) )^(m + 1) / (2 * pi * D^2);
    end
    
    H_LOS_W_t_all(i, :) = H_LOS_W_t;
    P_recep(i, :) = H_LOS_W_t * P_transmision;


    for j = 1:size(H_LOS_W_t, 2)
        Pr_ruido_termico_vec(i, j) = H_LOS_W_t(j) * P_transmision + thermal_noise;  % Potencia recibida al adicionar ruido térmico en el receptor
        Pr_ruido_disparo_vec(i, j) = H_LOS_W_t(j) * P_transmision + I_n;            % Potencia recibida al adicionar ruido de disparo en el receptor
        Pr_ruido_disparo_y_termico_vec(i, j) = H_LOS_W_t(j) * P_transmision + I_n + thermal_noise; 
        
        % Llamada a la función calculo_SNR para cada LED y almacenar el valor en SNR_Paralelos_por_LED
        P_recepcion(i,j) = H_LOS_W_t(j) * P_transmision; % Definir la potencia de recepción para la iteración
        SNR_Paralelos_por_LED(i, j) = calculo_SNR(k, T, B, R, q, I_b, P_recepcion(i,j));
        
    end
    
    fprintf('Iteración %d: Pr_ruido_termico=%.4e, %.4e, %.4e\n', i, Pr_ruido_termico_vec(i, 1), Pr_ruido_termico_vec(i, 2), Pr_ruido_termico_vec(i, 3));
    fprintf('Iteración %d: Pr_ruido_disparo_vec=%.4e, %.4e, %.4e\n', i, Pr_ruido_disparo_vec(i, 1), Pr_ruido_disparo_vec(i, 2), Pr_ruido_disparo_vec(i, 3));


    %Mostrar resultados de potencia al adicionar ruido térmico y de disparo
    fprintf('Iteración %d: Pr_ruido_disparo_y_termico_vec=%.4e, %.4e, %.4e\n', i, Pr_ruido_disparo_y_termico_vec(i, 1), Pr_ruido_disparo_y_termico_vec(i, 2), Pr_ruido_disparo_y_termico_vec(i, 3));

    Precepcion_t(i, :) = [h_i, d_cercana_i, H_LOS_W_t_all(i, 1), H_LOS_W_t_all(i, 2), H_LOS_W_t_all(i, 3), t_i];
    
    % Mostrar valores de Precepcion_t a medida que se llena
    fprintf('Iteración %d: h_i=%.2f, d_cercana_i=%.2f, H_los_1=%.4e, H_los_2=%.4e, H_los_3=%.4e, t_i=%.6e\n', ...
            i, Precepcion_t(i, 1), Precepcion_t(i, 2), Precepcion_t(i, 3), Precepcion_t(i, 4), Precepcion_t(i, 5), Precepcion_t(i, 6));


end

fprintf('\n\n');

disp('H_LOS_W_t de cada iteración:');
disp(H_LOS_W_t_all);

fprintf('\n       Valores de potencia cuando los Nt es paralelo a Nr \n');

disp('P_recep de cada iteración:');
disp(P_recep);


fprintf('\n       Valores de potencia y ruido cuando los Nt es paralelo a Nr \n');



disp('Pr_ruido_termico de cada iteración:');
disp((Pr_ruido_termico_vec));
fprintf('\n\n');
disp('Pr_ruido_disparo_vec de cada iteración:');
disp(Pr_ruido_disparo_vec);
fprintf('\n\n');
disp('Pr_ruido_disparo_y_termico_vec de cada iteración:');
disp(Pr_ruido_disparo_y_termico_vec);
fprintf('\n\n');


fprintf('------- Valores más altos de Pr_ruido_termico por iteración -------\n');
for i = 1:size(Pr_ruido_termico_vec, 1)
    max_value = max(Pr_ruido_termico_vec(i, :)); 
    fprintf('Iteración %d: Valor máximo=%.4e\n', i, max_value);
    min_value = min(Pr_ruido_termico_vec(i, :)); 
    fprintf('Iteración %d: Valor mínimo=%.4e\n\n', i, min_value);
end

fprintf('------- Valores más altos de Pr_ruido_disparo por iteración planos paralelos -------\n');
for i = 1:size(Pr_ruido_disparo_vec, 1)
    max_value = max(Pr_ruido_disparo_vec(i, :)); 
    fprintf('Iteración %d: Valor máximo=%.4e\n', i, max_value);
    min_value = min(Pr_ruido_disparo_vec(i, :)); 
    fprintf('Iteración %d: Valor mínimo=%.4e\n\n', i, min_value);
end


fprintf('------- Valores más altos de Pr_ruido_disparo_y_termico_vec por iteración planos paralelos -------\n');
for i = 1:size(Pr_ruido_disparo_y_termico_vec, 1)
    max_value = max(Pr_ruido_disparo_y_termico_vec(i, :)); 
    fprintf('Iteración %d: Valor máximo=%.4e\n', i, max_value);
    min_value = min(Pr_ruido_disparo_y_termico_vec(i, :)); 
    fprintf('Iteración %d: Valor mínimo=%.4e\n\n', i, min_value);
end


disp('----------------------------------------------------------------------------------------------------------');  %
disp('-----------------------Calculos de SNR al adicionar ruido termico y de disparo-----------------------');
disp('----------------------------------------------------------------------------------------------------------');  % 


disp('SNR de plano receptor paralelo al plano transmisor (SNR_Paralelos_por_LED):');
disp(SNR_Paralelos_por_LED);

x_values = 1:8; 

% Graficar los valores de SNR para los tres LEDs
figure; 
hold on; 
set(gcf, 'Color', 'w'); 

% Curva para el primer LED
plot(x_values, SNR_Paralelos_por_LED(:, 1), '-o', 'DisplayName', 'LED 1', 'LineWidth', 2);

% Curva para el segundo LED
plot(x_values, SNR_Paralelos_por_LED(:, 2), '-s', 'DisplayName', 'LED 2', 'LineWidth', 2);

% Curva para el tercer LED
plot(x_values, SNR_Paralelos_por_LED(:, 3), '-d', 'DisplayName', 'LED 3', 'LineWidth', 2);


xlabel('Posiciones en x', 'FontSize', 12); 
ylabel('SNR [dB]', 'FontSize', 12);       
title('Curvas de SNR por LED con Nt y Nr paralelos', 'FontSize', 14); 
legend('show'); 
grid on; 
hold off; 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




disp('----------------------------------------------------------------------------------');  %
disp('-------------Si suponemos que el FD NO está siempre paralelo al techo-------------');
disp('----------------------------------------------------------------------------------');  % 

disp('----Debemos calcular los valores de HLOS en función de los vectores de Nt y Nr----'); 


Psi = zeros(8, 3);                      % Almacenar valores de Psi para cada combinación i,j
distancia_euclidiana = zeros(8, 3);     % Almacenar distancias euclidianas
HLOS = zeros(8, 3);                     % Almacenar HLOS para cada combinación i,j
Pr_plano_no_paralelos = zeros(8, 3);    % Almacenar potencia por plano no paralelo

% Inicializar vectores de ruido
Pr_ruido_termico_vec_planos_no_paralelos = zeros(8, 3); 
Pr_ruido_disparo_vec_planos_no_paralelos = zeros(8, 3); 
Pr_ruido_disparo_y_termico_vec_planos_no_paralelos = zeros(8, 3); 


% Inicializaar vectores de SNR
SNR_NO_Paralelos_por_LED = zeros(8,3);
SNR_MAX_NO_Paralelos = zeros(8,2);      
SNR_MIN_NO_Paralelos = zeros(8,2);

% Inicialización del vector de potencia más alta por receptor
potencias_por_receptor = zeros(8, 1);
mejor_transmisor = zeros(8, 1);

for i = 1:8
    Xr = Pos_FD(i, 1);
    Yr = Pos_FD(i, 2);
    Zr = Pos_FD(i, 3);
    
    for j = 1:3
        Xt = LED_pos(j, 1);
        Yt = LED_pos(j, 2);
        Zt = LED_pos(j, 3);
        
        % Calcular Psi (ángulo entre el receptor i y el transmisor j)
        Psi(i, j) = acos((Xt * Xr + Yt * Yr + Zt * Zr) / ...
            (sqrt(Xr^2 + Yr^2 + Zr^2) * sqrt(Xt^2 + Yt^2 + Zt^2)));
        
        % Calcular distancia euclidiana entre el receptor i y el LED j
        distancia_euclidiana(i, j) = sqrt((Xt - Xr)^2 + (Yt - Yr)^2 + (Zt - Zr)^2);
        
        % Calcular HLOS para el receptor i y el LED j
        HLOS(i, j) = Ts * G_Con * Adet * (m + 1) * (cosd(theta)^m) * cos(Psi(i, j)) / ...
                     (2 * pi * distancia_euclidiana(i, j)^2);
        
        % Calcular potencia recibida por el plano no paralelo (Pr_plano_no_paralelos)
        Pr_plano_no_paralelos(i, j) = P_transmision * HLOS(i, j);
        
        % Calcular potencias con ruido térmico, de disparo y ambos
        Pr_ruido_termico_vec_planos_no_paralelos(i, j) = Pr_plano_no_paralelos(i, j) + thermal_noise;
        Pr_ruido_disparo_vec_planos_no_paralelos(i, j) = Pr_plano_no_paralelos(i, j) + I_n;
        Pr_ruido_disparo_y_termico_vec_planos_no_paralelos(i, j) = Pr_plano_no_paralelos(i, j) + I_n + thermal_noise;

        SNR_NO_Paralelos_por_LED(i, j) = calculo_SNR(k, T, B, R, q, I_b, Pr_plano_no_paralelos(i, j));

    end
    
    % Encontrar el LED con la mayor potencia para el receptor i
    [potencias_por_receptor(i), mejor_transmisor(i)] = max(Pr_plano_no_paralelos(i, :));
end

column_names = {'Receptor', 'Transmisor', 'Psi', 'Distancia Euclidiana', 'HLOS', 'Pr_plano_no_paralelos'};
result_data = [];

for i = 1:8
    for j = 1:3
        result_data = [result_data; i, j, Psi(i, j), distancia_euclidiana(i, j), HLOS(i, j), Pr_plano_no_paralelos(i, j)];
    end
end


result_table = array2table(result_data, 'VariableNames', column_names);
disp('Tabla de resultados (llegada_al_receptor):');
disp(result_table);

summary_data = [];
for i = 1:8
    summary_data = [summary_data; i, mejor_transmisor(i), Psi(i, mejor_transmisor(i)), distancia_euclidiana(i, mejor_transmisor(i)), HLOS(i, mejor_transmisor(i)), potencias_por_receptor(i)];
end

summary_table = array2table(summary_data, 'VariableNames', column_names);
disp('Tabla de Resumen - Potencia más alta por receptor:');
disp(summary_table);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

disp('Pr_plano_no_paralelos de cada iteración:');
disp((Pr_plano_no_paralelos));



disp('--------------------------------------------------------------------------------------------------------------');
disp('--- Resultados de potencia con ruido cuando el ángulo entre el receptor y cada uno de los LEDs no es ideal ---');
disp('--------------------------------------------------------------------------------------------------------------');

disp('Pr_ruido_termico_planos_no_paralelos de cada iteración:');
disp((Pr_ruido_termico_vec_planos_no_paralelos));

disp('Pr_ruido_disparo_vec_planos_no_paralelos de cada iteración:');
disp(Pr_ruido_disparo_vec_planos_no_paralelos);

disp('Pr_ruido_disparo_y_termico_vec_planos_no_paralelos de cada iteración:');
disp(Pr_ruido_disparo_y_termico_vec_planos_no_paralelos);


fprintf('------- Valores más altos de Pr_ruido_termico por iteración -------\n');
for i = 1:size(Pr_ruido_termico_vec_planos_no_paralelos, 1)
    max_value = max(Pr_ruido_termico_vec_planos_no_paralelos(i, :)); 
    fprintf('Iteración %d: Valor máximo=%.4e\n', i, max_value);
    min_value = min(Pr_ruido_termico_vec_planos_no_paralelos(i, :)); 
    fprintf('Iteración %d: Valor mínimo=%.4e\n\n', i, min_value);
end

fprintf('------- Valores más altos de Pr_ruido_disparo por iteración -------\n');
for i = 1:size(Pr_ruido_disparo_vec_planos_no_paralelos, 1)
    max_value = max(Pr_ruido_disparo_vec_planos_no_paralelos(i, :)); 
    fprintf('Iteración %d: Valor máximo=%.4e\n', i, max_value);
    min_value = min(Pr_ruido_disparo_vec_planos_no_paralelos(i, :)); 
    fprintf('Iteración %d: Valor mínimo=%.4e\n\n', i, min_value);
end

fprintf('---- Valores más altos de Pr_ruido_disparo_y_termico_vec por iteración cuando los planos Nt y Nr no son paraleros ----\n');
for i = 1:size(Pr_ruido_disparo_y_termico_vec_planos_no_paralelos, 1)
    max_value = max(Pr_ruido_disparo_y_termico_vec_planos_no_paralelos(i, :)); 
    fprintf('Iteración %d: Valor máximo=%.4e\n', i, max_value);
    min_value = min(Pr_ruido_disparo_y_termico_vec_planos_no_paralelos(i, :));
    fprintf('Iteración %d: Valor mínimo=%.4e\n\n', i, min_value);
end



disp('----------------------------------------------------------------------------------------------------------');  %
disp('-----------------------Calculos de SNR al adicionar ruido termico y de disparo-----------------------');
disp('----------------------------------------------------------------------------------------------------------');  % 

disp('SNR de plano receptor no paralelo al plano transmisor (SNR_NO_Paralelos_por_LED):');
disp(SNR_NO_Paralelos_por_LED);


% Graficar los valores de SNR para los tres LEDs
figure; 
hold on; 
set(gcf, 'Color', 'w'); 

% Curva para el primer LED
plot(x_values, SNR_NO_Paralelos_por_LED(:, 1), '-o', 'DisplayName', 'LED 1', 'LineWidth', 2);

% Curva para el segundo LED
plot(x_values, SNR_NO_Paralelos_por_LED(:, 2), '-s', 'DisplayName', 'LED 2', 'LineWidth', 2);

% Curva para el tercer LED
plot(x_values, SNR_NO_Paralelos_por_LED(:, 3), '-d', 'DisplayName', 'LED 3', 'LineWidth', 2);


xlabel('Posiciones en x', 'FontSize', 12); 
ylabel('SNR [dB]', 'FontSize', 12);       
title('Curvas de SNR por LED con Nt y Nr NO paralelos', 'FontSize', 14); 

legend('show'); 
grid on; 
hold off; 



disp('----------------------------------------------------------------------------------------------------------');  %
disp('------------Calculos de comparación SNR al juntar el mejor y peor caso en curvas de MAX y MIN-------------');
disp('----------------------------------------------------------------------------------------------------------');  % 



% Obtener los valores máximos y mínimos por fila (por cada iteración)
max_values_paralelos = max(SNR_Paralelos_por_LED, [], 2); % Máximos por fila (Paralelos)
min_values_paralelos = min(SNR_Paralelos_por_LED, [], 2); % Mínimos por fila (Paralelos)

max_values_no_paralelos = max(SNR_NO_Paralelos_por_LED, [], 2); % Máximos por fila (No Paralelos)
min_values_no_paralelos = min(SNR_NO_Paralelos_por_LED, [], 2); % Mínimos por fila (No Paralelos)


figure; 
hold on; 
set(gcf, 'Color', 'w'); 


% Graficar los valores máximos de los paralelos
plot(x_values, max_values_paralelos, '-o', 'DisplayName', 'Máximos planos paralelos', 'LineWidth', 2);

% Graficar los valores máximos de los no paralelos
plot(x_values, max_values_no_paralelos, '-x', 'DisplayName', 'Máximos planos No paralelos', 'LineWidth', 2);

xlabel('Posiciones en x', 'FontSize', 12); 
ylabel('SNR [dB]', 'FontSize', 12);       
title('Curva de SNR Máximos', 'FontSize', 14); 
legend('show'); 
grid on; 
hold off;

porcentaje_max_values = (abs (max_values_no_paralelos - max_values_paralelos) ./ max_values_paralelos) * 100;

disp('Porcentaje de max_values_paralelos sobre max_values_no_paralelos:');
for i = 1:length(porcentaje_max_values)
    fprintf('%.2f%%\n', porcentaje_max_values(i));
end


% Graficar el porcentaje de max_values_paralelos sobre max_values_no_paralelos
figure;
plot(x_values, porcentaje_max_values, '-o', 'DisplayName', 'Porcentaje Máximos', 'LineWidth', 2);
set(gcf, 'Color', 'w'); 

xlabel('Posiciones en x', 'FontSize', 12); 
ylabel('Porcentaje [%]', 'FontSize', 12);  
title('Porcentaje de Máximos Paralelos sobre No Paralelos', 'FontSize', 14); 
grid on; 




figure; 
hold on; 
set(gcf, 'Color', 'w'); 

% Graficar los valores mínimos de los paralelos
plot(x_values, min_values_paralelos, '-s', 'DisplayName', 'Mínimos planos paralelos', 'LineWidth', 2);

% Graficar los valores mínimos de los no paralelos
plot(x_values, min_values_no_paralelos, '-d', 'DisplayName', 'Mínimos planos No paralelos', 'LineWidth', 2);

xlabel('Posiciones en x', 'FontSize', 12); 
ylabel('SNR [dB]', 'FontSize', 12);       
title('Curva de SNR Mínimos', 'FontSize', 14); 
legend('show');
grid on;  
hold off; 



porcentaje_min_values = (abs (min_values_no_paralelos - min_values_paralelos) ./ min_values_paralelos) * 100;

disp('Porcentaje de min_values_paralelos sobre min_values_no_paralelos:');
for i = 1:length(porcentaje_min_values)
    fprintf('%.2f%%\n', porcentaje_min_values(i));
end

% Graficar el porcentaje de min_values_paralelos sobre min_values_no_paralelos
figure;
set(gcf, 'Color', 'w'); 
plot(x_values, porcentaje_min_values, '-o', 'DisplayName', 'Porcentaje Máximos', 'LineWidth', 2);

xlabel('Posiciones en x', 'FontSize', 12); 
ylabel('Porcentaje [%]', 'FontSize', 12);  
title('Porcentaje de mínimos Paralelos sobre No Paralelos', 'FontSize', 14); % Título de la gráfica
grid on; 



                        %%% ---------- DEFINICIÓN DE FUNCIONES ---------- %%%%

                        
                        
function SNR = calculo_SNR(k, T, B, R, q, I_b, P_recepcion)

    % Cálculo del ruido térmico en potencia (watts)
    V_rms_thermal = sqrt(4 * k * T * R * B);  % Voltaje RMS del ruido térmico
    ruido_termico = (V_rms_thermal^2) / R;    % Potencia del ruido térmico en watts

    % Alternativamente: ruido_termico = (4 * k * T * B) / R;
    
    % Cálculo del ruido de disparo en potencia (watts)
    I_rms_shot = sqrt(2 * q * I_b * B);  % Corriente RMS del ruido de disparo
    ruido_disparo = (I_rms_shot^2) * R;  % Potencia del ruido de disparo en watts
    
    % Alternativamente: ruido_disparo = 2 * q * I_b * B * R;
    
    % Cálculo del ruido total en potencia (watts)
    ruido_total = ruido_disparo + ruido_termico;
    
    % Cálculo del SNR (relación señal/ruido en términos de potencia)
    SNR = 10 * log10( ( P_recepcion) / ruido_total);  % relación en términos de potencia

end



disp('---------------------------------------------------------------------------------');  %
disp('------------Calculos de distancia de error de localización en metros-------------');
disp('---------------------------------------------------------------------------------');  % 

% Inicialización: convertir matriz de distancias a 2D
dist_estimada = matriz_dist_ecuclidiana'; % Cada fila contiene las distancias estimadas

% Inicialización para guardar las soluciones (x, y)
pos_estimada = zeros(size(dist_estimada, 1), 2);
error_distancia = zeros(8, 1);

%pos_estimada = calcular_posiciones(dist_estimada, LED_pos);
pos_estimada = calcularPosiciones(dist_estimada, LED_pos);
disp('Posiciones estimadas en casos sin ruido:');
disp(pos_estimada);


% Coordenadas reales (Pos_real) en 2D
Pos_real = zeros(8, 2);
for i = 1:8 % En total son 8 puntos
    Pos_real(i, :) = [i, 2]; % (x, y)
end


error_distancia = calcular_error(Pos_real, pos_estimada);
disp('Error de distancia en casos sin ruidos:');
disp(error_distancia);

plotPositions(Pos_real, pos_estimada);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_i = 3.6; % Altura del receptor
z_t = 2.6; % Altura del transmisor

d_estimada_no_paralelos = zeros(8, 3);
Pos_est_no_paralelos = zeros(8,2);
error_no_paralelos = zeros(8,1);

% Cálculo de d_estimada_no_paralelos
for i = 1:size(Pr_ruido_disparo_y_termico_vec_planos_no_paralelos, 1)
    for j = 1:size(Pr_ruido_disparo_y_termico_vec_planos_no_paralelos, 2)

        P_LOS_i = Pr_ruido_disparo_y_termico_vec_planos_no_paralelos(i, j);
        d_estimada_no_paralelos(i, j) = ((Ap * (m + 1) / (2 * pi)) * fov * G_Con * P_transmision * ((z_i - z_t)^(m + 1) / P_LOS_i))^(1 / (m + 3));
   
    end
end


disp('Matriz d_estimada_no_paralelos:');
disp(d_estimada_no_paralelos);

Pos_est_no_paralelos = calcularPosiciones(d_estimada_no_paralelos, LED_pos);
disp('Posiciones estimadas (No Paralelos):');
disp(Pos_est_no_paralelos);

error_no_paralelos = calcular_error(Pos_real, Pos_est_no_paralelos);
disp('Error de distancia (No Paralelos):');
disp(error_no_paralelos);

plotPositions(Pos_real, Pos_est_no_paralelos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_estimada_paralelos = zeros(8, 3);
Pos_est_paralelos = zeros(8,2);
error_paralelos = zeros(8,1);

% Cálculo de d_estimada_paralelos
for i = 1:size(Pr_ruido_disparo_y_termico_vec, 1)
    for j = 1:size(Pr_ruido_disparo_y_termico_vec, 2)

        P_LOS_i = Pr_ruido_disparo_y_termico_vec(i, j);
        d_estimada_paralelos(i, j) = ((Ap * (m + 1) / (2 * pi)) * fov * G_Con * P_transmision * ((z_i - z_t)^(m + 1) / P_LOS_i))^(1 / (m + 3));
    
    end
end

disp('Matriz d_estimada_paralelos:');
disp(d_estimada_paralelos);

Pos_est_paralelos = calcularPosiciones(d_estimada_paralelos, LED_pos);
disp('Valores estimados de X e Y para cada iteración:');
disp(Pos_est_paralelos);


error_paralelos = calcular_error(Pos_real, Pos_est_paralelos);
disp('Error de distancia (Paralelos):');
disp(error_paralelos);

plotPositions(Pos_real, Pos_est_paralelos);

% Crear una figura para comprender como se aplica TOA
figure;
hold on; 
c = 3e8; 
tiempos_llegada = d_estimada_paralelos / c;
set(gcf, 'Color', 'w'); 
colores = ['r', 'g', 'b']; % Rojo para LED 1, Verde para LED 2, Azul para LED 3

for i = 1:size(LED_pos, 1)
    line([tiempos_llegada(1, i) tiempos_llegada(1, i)], [0 d_estimada_paralelos(1, i)], 'Color', colores(i), 'LineWidth', 2);
end

% Añadir etiquetas y leyenda
xlabel('Tiempo de Llegada (s)');
ylabel('Distancia Estimada (metros)');
title('Representación de Impulsos de Llegada para cada LED');
legend({'LED 1', 'LED 2', 'LED 3'}, 'Location', 'northeast');
grid on;
hold off; % Liberar el gráfico

%%%%%%%%%%%%%%%%%%%%%% FUNCIONES DE LOCALIZACION %%%%%%%%%%%%%%%%%%%%%%
function Pos_est = calcularPosiciones(d_estimada, LED_pos)

    Pos_est = zeros(size(d_estimada, 1), 2); % N filas, 2 columnas para [X, Y]

    % Iterar sobre cada fila de d_estimada para calcular X e Y
    for i = 1:size(d_estimada, 1)
        XLed = LED_pos(:, 1);
        YLed = LED_pos(:, 2);

        % Extraer las distancias estimadas para la iteración actual
        dist_actual = d_estimada(i, :);

        % Ordenar las distancias y seleccionar las dos menores
        [~, indices_ordenados] = sort(dist_actual);
        idx1 = indices_ordenados(1); % Menor distancia
        idx2 = indices_ordenados(2); % Segunda menor distancia

        % Definir distancias y coordenadas asociadas
        dest1 = dist_actual(idx1);
        dest2 = dist_actual(idx2);
        Xled1 = XLed(idx1);
        Xled2 = XLed(idx2);

        % Calcular X usando la fórmula dada
        numerator = dest1^2 - dest2^2 - Xled1^2 + Xled2^2;
        denominator = -2 * Xled1 + 2 * Xled2;

        if denominator ~= 0
            X = abs(numerator / denominator);
        else
            X = NaN; % Si hay división por cero, asignar NaN
        end

        % Fijar Y como constante (por simplificación)
        Y = encontrar_y_optimo(X, Xled1, dest1);

        % Guardar los valores en la matriz de posiciones estimadas
        Pos_est(i, :) = [X, Y];
    end
end

function y_opt = encontrar_y_optimo(x, x_led, dest)
    % Función para encontrar el valor óptimo de y en el rango [1, 3]
    % que minimiza la ecuación (x - x_led)^2 + (y - 2)^2 - dest^2

    % Definir la función objetivo
    objective = @(y) abs((x - x_led)^2 + (y - 4)^2 - dest^2);

    % Configurar opciones para fminbnd con mayor precisión
    options = optimset('TolX', 23e-6, 'Display', 'off');

    % Usar fminbnd para encontrar el valor óptimo de y en el rango [1, 3]
    y_opt = fminbnd(objective, 1.25, 2.75, options);
end
%%%%%%%%%%%%%%%%%%%%% Fin de calculo de localización %%%%%%%%%%%%%%%%%%%%%




function error_distancia = calcular_error(Pos_real, pos_estimada)
    error_distancia = zeros(size(Pos_real, 1), 1);

    % Calcular la distancia para cada punto
    for i = 1:size(Pos_real, 1)
        % Extraer coordenadas reales (x, y) y estimadas (x, y)
        x_real = Pos_real(i, 1);
        y_real = Pos_real(i, 2);

        x_est = pos_estimada(i, 1);
        y_est = pos_estimada(i, 2);

        % Calcular la distancia euclidiana
        error_distancia(i) = sqrt((x_real - x_est)^2 + (y_real - y_est)^2);
    end
end



function plotPositions(real_pos, estim_pos)

    % Verificar dimensiones de entrada
    if size(real_pos, 2) ~= 2 || size(estim_pos, 2) ~= 2
        error('Ambos vectores deben ser matrices Nx2 con columnas [x, y].');
    end


    figure;
    hold on;
    set(gcf, 'Color', 'w'); 
    plot(real_pos(:, 1), real_pos(:, 2), 'kx', 'DisplayName', 'Posición real'); % Cruces negras
    plot(estim_pos(:, 1), estim_pos(:, 2), 'ro', 'DisplayName', 'Posición estimada'); % Círculos rojos
    legend;
    xlabel('X');
    ylabel('Y');
    axis([0 8 0 4]); % Ajustar el rango del gráfico
    axis equal;
    grid on;
    title('Desviación posición real v/s estimada');
    hold off;
end
