function [h, spyr] = spyrDisp4(x, pind, P, varargin)
    
    S = P(1);
    O = P(2);
    N = P(3);
    
    I = zeros(S*N, (O+1)*N);
    
    if nargin >= 4;
        smoother = varargin{1};
    else
        smoother = 0;
    end
    if nargin >= 5;
        sd = varargin{2};
        
        zf = 2.^(log2(N/4)-(S:-1:0))/sd;
        for s = 1:S
            zfl{s} = sprintf('%2.1f-%2.1f', zf(s), zf(s+1));
        end
    else
        for s = 1:S
            zfl{s} = sprintf('%i', s);
        end
    end

    z = 0;
    for r = 1:size(pind, 1);
        zNext = prod(pind(r,:));

        f = x(z+1:z+zNext);
        f = imresize(reshape(f, pind(r,:)), (N)*[1 1], 'nearest');
        spyr(1+(r-1)*N^2:(r)*N^2) = f;        

        z = z+zNext;
    end
    x = spyr;

    %isolate info
    fhr = x(1:N^2);
    for s = 1:S
        fb{s} = x(((s-1)*O+1)*N^2+1:((s*O+1)*N^2));
    end
    flr = x((S*O+1)*N^2+1:(S*O+2)*N^2);

    I = [];
    for s = 1:S
        for o = 1:O

            f = reshape(fb{s}((o-1)*N^2+1:o*N^2), [N N]);

            if smoother ~=0
               dx = -32:1:32;
               fx = 1/(sqrt(2*pi*smoother^2)) .* exp(-1*dx.^2/(2*(smoother^2)));
               fx2 = fx'*fx;
               f = imfilter(f, fx2);
            end
            I((S-s)*N+1:(S-s+1)*N, (o-1)*N+1:o*N) = f;
        end
    end

    %append residuals
    I(1:N, O*N+1:(O+1)*N) = reshape(flr, [N N]);
    I(end-N+1:end, O*N+1:(O+1)*N) = reshape(fhr, [N N]);

    I = I / 2 + 0.5;
    I(I>1) = 1;
    I(I<0) = 0;

    for s = 1:S-1
        I(s*N:s*N+1,1:end) = 1;
    end
        
    oriLabelSet = mod(90:-180/O:-90+1, 180);
    oriLabel = {};
    for o = 1:O
        I(1:end,N*o:N*o+1) = 1;
        oriLabel{o} = sprintf('%3i', oriLabelSet(o));
        if oriLabelSet(o) == 90
            oriLabel{o} = '90 (vert.)';
        elseif oriLabelSet(o) == 0
            oriLabel{o} = '0 (horiz.)';
        end
    end
    oriLabel{O+1} = 'Res.';
    
    scaleLabel = zfl;
    scaleLabel{1} = ['Coarse ' num2str(scaleLabel{1})];
    scaleLabel{end} = ['Fine ' num2str(scaleLabel{end})];

    
    h = imagesc(1:(O+1)*N, 1:(S)*N, I);
    axis equal tight;
    set(gca, 'tickdir', 'out', 'xtick', N/2:N:(O+1)*N, 'xticklabel', oriLabel, 'ytick', N/2:N:S*N, 'yTickLabel', scaleLabel);
    colormap(gray(256))
    %turn I back into a vector and outputz
end

