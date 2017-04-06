function expm_stats()

% Statistics to compute expm.
% --------------------------
% procedure to evaluate the mean computational time of expm, for 
% random matrices generated with increasing size or sigma.

num_samples_per_steps_size = 50;
list_size = 200; 
time_increasing_size = zeros(list_size, num_samples_per_steps_size);


num_samples_per_steps_sigma = 50;
list_sigma = 200;
fixed_size_for_sigma = 5;
time_increasing_sigma = zeros(list_sigma, num_samples_per_steps_sigma);


for siz=1:list_size
   for sample=1:num_samples_per_steps_size

       A = randn(siz,siz);
       tic
       expm(A);
       time_increasing_size(siz, sample) = toc;

   end
end


for sigma=1:list_sigma
   for sample=1:num_samples_per_steps_sigma

       A = sigma * randn(fixed_size_for_sigma, fixed_size_for_sigma);
       tic
       expm(A);
       time_increasing_sigma(sigma, sample) = toc;

   end
end

g = figure();

set(g, 'Position', [100, 100, 1150, 850]);

figure(g)

subplot(211)

x = 1:size(time_increasing_size,1) ;
y = mean(time_increasing_size, 2)';
y_std = std(time_increasing_size, 1, 2);

hold on
H = shadedErrorBar(x, y, y_std, '-g', 0);
hold off

legend([H.mainLine, H.patch], ...
     '\mu', '\sigma', ...
     'Location', 'Northwest');

title('time to compute expm increasing size')
xlabel('size of the squared matrix')
ylabel('time to compute expm (sec.)')

subplot(212)
x = 1:size(time_increasing_sigma,1);
y = mean(time_increasing_sigma, 2)';
y_std = std(time_increasing_sigma, 1, 2);


hold on
H = shadedErrorBar(x, y, y_std, '-g', 0);
hold off

legend([H.mainLine, H.patch], ...
     '\mu', '\sigma', ...
     'Location', 'Northwest');

title('time to compute expm increasing sigma')
xlabel('sigma of the random squared matrix')
ylabel('time to compute expm (sec.)')

end


