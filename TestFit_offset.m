
function cost = TestFit_offset( params, patch_sig, xsens_sig, t_expand )
    t_offset = params(1);
    
    new_patch_t = ( patch_sig(:,1) ) * t_expand + t_offset;
    
    xsens_comp = interp1( xsens_sig(:,1), xsens_sig(:,2), new_patch_t, 'spline' );
    
    cost = sum( ( xsens_comp - patch_sig(:,2) ) .* ( xsens_comp - patch_sig(:,2) ) );
end