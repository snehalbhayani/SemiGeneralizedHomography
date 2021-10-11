
function d = smoothData(d, wsiz)

    if wsiz==1
        
        return;
    end

    d = [d(1)*ones(1, wsiz) d d(end)*ones(1, wsiz)];
    
	ker = fspecial('gaussian',[1 wsiz], 4);

    b = conv(d, ker);
	d = b((wsiz+floor(wsiz/2)):(floor(end-wsiz/2)-wsiz));
end

