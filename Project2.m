%% Define a triangular pulse shaping function and Random bit stream generation
triangular_pulse = [5 4 3 2 1] / sqrt(55); % triangluar pulse to be convolved
                                           % with the bitstream.

% The pulse is normalized so that its energy is = 1,ensuring a fair 
% comparison of filter performance% 

bit_stream = randi([0, 1], 1, 10); %creates Random array of 10 bits, 
                                   % The input data to be transmitted.

Polar_NRZ_stream = 2*bit_stream - 1;   %convert bitstream into polar NRZ line code.
%% Upsampling the bit stream
upsampled_stream = upsample(Polar_NRZ_stream, 5); %Inserts 4 zeros between each bit, 
                                             % resulting in a sampling rate of
                                             % 5 samples per symbol.

% Ensures that each bit lasts 1 second, and sampling occurs every 200 ms.%
%% Define Time Reference for Bit Stream
signal_length = length(upsampled_stream); % Calculate the length of 
                                         % the upsampled signal.

% Define a time axis where each sample is equally spaced 200 ms apart.                                         
stream_time_reference = linspace(0, (signal_length - 1) * 0.2, signal_length);
%% Pulse Shaping (T.D convolution <-> F.D Multiplication)
convolved_stream = conv(upsampled_stream,triangular_pulse);
%% Define Time Reference for Convolved Signal
convolved_stream_length = length(convolved_stream);
convolved_time_reference = linspace(0, (convolved_stream_length - 1) * 0.2, convolved_stream_length);

%% Plot the Bit Stream, Pulse Shaping Function and the Convolved Signal
figure;
stem(stream_time_reference, upsampled_stream); % stem plot to visualize the upsampled bit stream.
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Bit Stream Samples');

figure;
plot(triangular_pulse);
xlabel('Sample Index'); 
ylabel('Amplitude'); 
title('Triangular Pulse ');

figure;
stem(convolved_time_reference, convolved_stream); % stem plot to visualize the upsampled bit stream.
xlabel('Time (seconds)'); 
ylabel('Amplitude'); 
title('Convolved Output Samples');

%% Matched Filter implementation, as it maximizes the SNR at sampling instants.
matched_filter = fliplr(triangular_pulse); % generating the matched filter,
                                         % by time-reversing the pulse
matchedfilter_output = conv(convolved_stream, matched_filter);
matchedfilter_time_reference = linspace(0, (length(matchedfilter_output) - 1) * 0.2, length(matchedfilter_output));

%% Define Rect function and normalize it
rect = [1 1 1 1 1]/sqrt(5); % Normalizes rect pulse filter
rect_output = conv(convolved_stream, rect); %Convolve the received signal with the rect filter.
rect_time_reference = linspace(0, (length(rect_output) - 1) * 0.2, length(rect_output));

%% Requirement 1(a): Plot Matched Filter vs. Rectangular Filter 
figure;
subplot(2,1,1);
plot(matchedfilter_time_reference, matchedfilter_output);
hold on;
stem(matchedfilter_time_reference, matchedfilter_output, 'green');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Output of the Matched Filter');

subplot(2,1,2);
plot(rect_time_reference, rect_output, 'red');
hold on;
stem(rect_time_reference, rect_output, 'red');
xlabel('Time (seconds)');
ylabel('Amplitude'); 
title('Output of the Rectangular Pulse');



%% Correlator Implementation
% Preallocate output
correlator_output = zeros(1, length(convolved_stream)); 

% Iterate over symbols (5 samples per symbol)
for i = 1:5:length(convolved_stream)-4  
    sum = 0;
    for j = 0:4 
        % Compute dot product (correlation)
        sum = sum + convolved_stream(i+j) * triangular_pulse(j+1);
    end
    % Store the final accumulated sum at the last sample of the symbol
    correlator_output(i+4) = sum; 
end


correlator_time_reference = linspace(0, (length(correlator_output) - 1) * 0.2, length(correlator_output));
%% Plot Correlator vs. Matched Filter Output
figure;
subplot(2,1,1);
plot(correlator_time_reference, correlator_output, 'red');
hold on;
stem(correlator_time_reference, correlator_output, 'red');
hold on;
legend('Output of the correlator', 'Sampled output of the correlator');
xlabel('Time (seconds)'); 
ylabel('Amplitude'); 
title("Correlator's Output");

% Matched filter
subplot(2,1,2);
plot(matchedfilter_time_reference, matchedfilter_output, 'green');
hold on;
stem(matchedfilter_time_reference, matchedfilter_output, 'green');
xlabel('Time (seconds)'); 
ylabel('Amplitude'); 
legend('Output of the matched filter', 'Sampled output of the matched filter');
title("Matched Filter's Output");

%% Requirement 1(b): Plot Correlator vs. Matched Filter on the Same Plot
figure;
plot(matchedfilter_time_reference, matchedfilter_output, 'b'); % Matched filter in blue
hold on;
plot(correlator_time_reference, correlator_output, 'r'); % Correlator output in red
stem(matchedfilter_time_reference(1:5:end), matchedfilter_output(1:5:end), 'bo'); % Sampled points for Matched filter
stem(correlator_time_reference(1:5:end), correlator_output(1:5:end), 'ro'); % Sampled points for Correlator
xlabel('Time (seconds)'); 
ylabel('Amplitude'); 
legend('Matched Filter Output', 'Correlator Output', 'Sampled Matched Filter Output', 'Sampled Correlator Output');
title("Comparison of Matched Filter and Correlator Outputs");
