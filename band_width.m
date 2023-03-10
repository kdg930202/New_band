function [result,center] = band_width(locate_bloch_22,omega,omegal)



band_gap = find(locate_bloch_22(2,:)==0);



for j = 1:length(band_gap)
%     band_0_list(1) = band_gap(2);

    if band_gap(1) ~= 1
        
        if band_gap(j+1) - band_gap(j) == 1
    %         display(i)
           band_0_list(j) = band_gap(j); 
    %         band_0_list(i) = band_gap(i+1);
        else
           band_0_list(j) = band_gap(j);
           break
        end
    else
        band_0_list(1) = band_gap(2);
        if band_gap(j+2) - band_gap(j+1) == 1
            band_0_list(j+1) = band_gap(j+2);     
        else
%             band_0_list(i+1) = band_gap(i);
            break 
        end
    end
end


ome_list = omega/omegal;
result = ome_list(band_0_list(end)) - ome_list(band_0_list(1));
center = (ome_list(band_0_list(end)) + ome_list(band_0_list(1)))/2;


end