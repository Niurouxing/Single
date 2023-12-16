clear
clc



SampleNum   = 100000;       % max sample no.
iterNum     = 30;% 7;        % iter no.   
delta       = 0.3;          % damping factor
delta2      = 0.3;
rho         = 0.3;
max_nberr   = 1000;

TxAntNum    = 16;
RxAntNum    = 32;
ModType 	= 4; 		
isCorr = 0;
% norm_f      = 1 / sqrt(10);  % normalize the Transimite Symnbol 

switch(ModType)
    case 2
        norm_f=1/sqrt(2);
        QAM=4;
        s=[ -1 +1 ]';
    case 4
        norm_f = 1/sqrt(10);
        QAM=16;
        s=[ -3 -1 +1 +3  ]';
    case 6
        norm_f=1/sqrt(42);
        QAM=64;
        s=[-7  -5 -3 -1 +1 +3 +5 +7 ]';
    case 8
        QAM=256;
        norm_f=1/sqrt(170);
        s=[-15 -13 -11 -9 -7  -5 -3 -1 +1 +3 +5 +7 +9 +11 +13 +15]';
end
N 			= 1; 						% Block length
InfoLen 	= N * ModType * TxAntNum; 	% Information length

EbN0db      = 6:2:16;
SampleTest  = zeros(length(EbN0db),1);  

H           = zeros(2*RxAntNum,2*TxAntNum);  
RxSymbol    = zeros(2*RxAntNum,1); % Received
TxSym_Detected = zeros(2*TxAntNum, 1);

slen = length(s);
sym = s*norm_f; % normalized symbol set

error_real = zeros(iterNum, length(EbN0db));
error_sym = zeros(iterNum, length(EbN0db));
error_MMSE = zeros(1, length(EbN0db));
error_s_MMSE = zeros(1, length(EbN0db));

for nEN =1:length(EbN0db)
    
    % initialization 
    SNR = EbN0db(nEN)
    
    loop = 0;
%     while loop == 0
    while error_real(iterNum, nEN) < max_nberr && loop < SampleNum
%     while error_real(iterNum, nEN) < max_nberr
        loop = loop + 1;
        % Data send via the channels
%         TxBits = zeros(InfoLen,1);
        TxBits = randi([0, 1], [InfoLen, 1]);
        Tx_mat = reshape(TxBits,ModType,InfoLen/ModType);
        Tx_H   = bi2de(Tx_mat.','left-msb');
        TxSymbol = GrayMapQAM((TxBits).', ModType).';
        
        % Channel Matrix
        if (isCorr)
            HMat = channel_mimo(RxAntNum,TxAntNum,rho); %Correlated channel
        else
		    HMat = sqrt(0.5) *(randn(RxAntNum, TxAntNum) + sqrt(-1)*randn(RxAntNum, TxAntNum)); %Rayleigh channel
        end
       
        % Complex channel matrix to real
        H = [real(HMat), -imag(HMat); imag(HMat), real(HMat)];
        
        % Received and sigma^2
        [RxSymbol_c, Nv] = acwgn_EbN0(HMat*TxSymbol, 1, SNR, TxAntNum, RxAntNum, ModType, 1);	% Add Gaussian noise
        Nv_MMSE = Nv;
        RxSymbol(1:RxAntNum,1)=real(RxSymbol_c);
        RxSymbol(RxAntNum+1:end,1)=imag(RxSymbol_c);

        % Complex to real for transmitted signal
        TxSymbol_r_0 = zeros(2*TxAntNum, 1);
        TxSymbol_r_0(1:TxAntNum)=1/norm_f*real(TxSymbol);
        TxSymbol_r_0(TxAntNum+1:end)=1/norm_f*imag(TxSymbol);
        Tx_B = symboltobits(TxSymbol_r_0, 2*TxAntNum, QAM); 
        TxSymbol_r_0 = round(TxSymbol_r_0);


        % Mutually Independent Q(xi) -- Parameter
        Alpha = 0.001*ones(2*TxAntNum, 1);
        sigma2_p  = zeros(2*TxAntNum, 1);
        mu_p = zeros(2*TxAntNum, 1);
	
        ytilde = H'*RxSymbol/Nv;  A = H'*H/Nv;
        Sigma_q = ( A+diag(Alpha)) \ eye(2*TxAntNum);
        Mu_q = Sigma_q * ytilde;
	    TxSym_Det = HD_EP_all(Mu_q/norm_f,s);
        Rx_B_MMSE = symboltobits(TxSym_Det, 2*TxAntNum,QAM);
        error_MMSE(nEN) = error_MMSE(nEN)+nnz(Rx_B_MMSE-Tx_B);
        Rx_S_MMSE = TxSym_Det(1:TxAntNum)*norm_f+TxSym_Det(TxAntNum+1:end)*norm_f*1j;
        error_s_MMSE(nEN) = error_s_MMSE(nEN)+nnz(TxSymbol-Rx_S_MMSE);
	    
       
           
        %Q(xi) -- Statistics
        sigma2 = 1*diag(Sigma_q);
        mu = Mu_q;
        
        % Q(xi)\i -- Statistics
        h2 = sigma2./(1-sigma2.*Alpha);
        t  = h2 .* (mu./sigma2); 
        
        % Convergence V
%         temp = (sum(A.^2,2)-diag(A).^2)./(diag(A).^3);
        V = diag(A);        
        
        %Initialization        
        ttt_new = t;   
        h2_new = h2;
        % Updated Q(xi)\i -- Statistics    
        ttt = t;         
     
        for k = 1:iterNum
            ttt = delta * ttt_new + (1-delta) * ttt;  
            h2  = delta2 * h2_new  + (1-delta2) * h2;
            % Hard Decision
%             for i = 1:2*TxAntNum 
%                 mu_p(i) = HD_EP_all(ttt(i)/norm_f,s); 
%                 mu_p(i) = mu_p(i)*norm_f; 
%             end

            value = -0.5*abs(ones(2*TxAntNum,1)*sym' - ttt*ones(1,sqrt(QAM))).^2./(h2*ones(1,sqrt(QAM)));
            dex = max(value,[],2);
            prob_muti = exp(value - dex*ones(1,sqrt(QAM)));  %calculate the probabilities used h and t
            prob_muti = prob_muti./(sum(prob_muti,2)*ones(1,sqrt(QAM)));  % normalization
            
            mu_p = prob_muti * sym; 
            var_temp = sum(prob_muti .* abs(ones(2*TxAntNum,1)*sym' - mu_p*ones(1,sqrt(QAM))).^2,2);
            for i = 1:2*TxAntNum
                var_temp(i) = max( [var_temp(i) 5e-7]);
            end
       
            m = ytilde-A*mu_p;
            ttt_new = m./V + mu_p;
            h2_new = var_temp;
%             h2_new = h2;
            
            %Cumulative Real-Time Bit Error
            TxSym_Detected = HD_EP_all(mu_p/norm_f,s);
            Rx_B_bp = symboltobits(TxSym_Detected, 2*TxAntNum,QAM);
            Rx_S_bp = TxSym_Detected(1:TxAntNum)*norm_f+TxSym_Detected(TxAntNum+1:end)*norm_f*1j;
            error_real(k, nEN) = error_real(k, nEN) + nnz(Rx_B_bp-Tx_B);
            error_sym(k,nEN) = error_sym(k,nEN)+nnz(Rx_S_bp-TxSymbol);
            
        end  

    end
    SampleTest(nEN) = loop; %SampleNum;
end

error_sym = error_sym ./ repmat((TxAntNum .* SampleTest'),iterNum,1);
error_real = error_real ./ repmat((TxAntNum * ModType.* SampleTest'),iterNum,1);
error_MMSE = error_MMSE ./ (TxAntNum * ModType.* SampleTest');
error_s_MMSE = error_s_MMSE ./ (TxAntNum .* SampleTest');

figure;
semilogy(EbN0db, error_MMSE, 'ko-',EbN0db, error_real(iterNum, :), 'r>-','LineWidth',2);
grid on;
legend('MMSE','FC-GAI-BP, iter=40');
xlabel('Average Received SNR(dB)')
ylabel('BER')