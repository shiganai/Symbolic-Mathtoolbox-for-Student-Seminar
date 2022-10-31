function [Coeffs, Terms] = coeffs_Vector(sym_Vector,target_Variables)

if ~iscolumn(sym_Vector)
    error('sym_Vector must be a column vector')
end

target_Variables = [target_Variables, 1];

Coeffs = sym(zeros(size(sym_Vector, 1), size(target_Variables, 2)));
Terms = target_Variables;

for jj = 1:size(sym_Vector, 1)
    [coeffs_Tmp, terms_Tmp] = coeffs(sym_Vector(jj), target_Variables(1:end-1));

    if jj == 1
        Coeffs = coeffs_Tmp;
        Terms = terms_Tmp;
    else

        coeffs_adding = sym( zeros( 1, size(Coeffs,2) ) );
        for ii = 1:size(terms_Tmp, 2)
            is_contained = false;

            for kk = 1:size(Terms,2)
                if isequal( terms_Tmp(ii), Terms(kk) )
                    % When terms_Tmp(ii) already exist in Terms
                    % add coeffs
                    coeffs_adding(kk) = coeffs_Tmp(ii);
                    is_contained = true;
                    break
                end
            end

            if ~is_contained
                % When terms_Tmp(ii) does not exist in Terms
                % add zero column at the end of Coeffs, and terms_Tmp(ii) at the end of Terms
                Coeffs = [Coeffs, zeros( size(Coeffs, 1) )];
                Terms = [Terms, terms_Tmp(ii)];
                coeffs_adding = [coeffs_adding, coeffs_Tmp(ii)];
            end
        end
        Coeffs = [Coeffs; coeffs_adding];
    end
end

end