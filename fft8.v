`timescale 1ns/1ps
module fft8(
    input  wire             clk,
    input  wire             rst,
    input  wire             valid_in,
    input  wire signed [63:0] xr_flat,
    input  wire signed [63:0] xi_flat,
    output reg              valid_out,
    output reg signed [63:0] Yr_flat,
    output reg signed [63:0] Yi_flat
);

    integer n;
    localparam POINTS = 8;
    localparam SCALE = 5;
	 
    reg signed [7:0] xr_arr [0:POINTS-1];
    reg signed [7:0] xi_arr [0:POINTS-1];
    always @(*) begin
        for(n=0;n<POINTS;n=n+1) begin
            xr_arr[n] = xr_flat[8*(n+1)-1 -: 8];
            xi_arr[n] = xi_flat[8*(n+1)-1 -: 8];
        end
    end

    // Twiddle factors for 8-point FFT
    reg signed [7:0] W_real [0:3];
    reg signed [7:0] W_imag [0:3];
    initial begin
        W_real[0]=32; W_imag[0]=0;    
        W_real[1]=0;  W_imag[1]=-32;  
        W_real[2]=23; W_imag[2]=-23;  
        W_real[3]=-23; W_imag[3]=-23; 
    end

    reg signed [7:0] stage1_r [0:POINTS-1], stage1_i [0:POINTS-1];
    reg signed [7:0] stage2_r [0:POINTS-1], stage2_i [0:POINTS-1];
    reg signed [7:0] stage3_r [0:POINTS-1], stage3_i [0:POINTS-1];
    reg valid1_r, valid2_r, valid3_r;
    reg valid1_i, valid2_i, valid3_i;

    reg signed [15:0] temp_r, temp_i;
    reg signed [7:0] mul_r, mul_i;
    reg [1:0] twiddle_idx;
	 
    integer a,b,c,d;
    reg signed [15:0] diff_real, diff_imag;
    reg signed [7:0] rot_r, rot_i;

	 //Stage 1 : 2-point butterflies
    always @(posedge clk or posedge rst) begin
        if(rst) begin
            valid1_r <= 0;
            for(n=0;n<POINTS;n=n+1) stage1_r[n]<=0;
        end else if(valid_in) begin
            valid1_r <= valid_in;
            for(n=0;n<POINTS/2;n=n+1) begin
                stage1_r[n]      <= xr_arr[n] + xr_arr[n+4];
                stage1_r[n+4]    <= xr_arr[n] - xr_arr[n+4];
            end
        end
    end
	 
    always @(posedge clk or posedge rst) begin
        if(rst) begin
            valid1_i <= 0;
            for(n=0;n<POINTS;n=n+1) stage1_i[n]<=0;
        end else if(valid_in) begin
            valid1_i <= valid_in;
            for(n=0;n<POINTS/2;n=n+1) begin
                stage1_i[n]      <= xi_arr[n] + xi_arr[n+4];
                stage1_i[n+4]    <= xi_arr[n] - xi_arr[n+4];
            end
        end
    end
	 
	 // Stage 2 : 4-point butterflies
	 always @(posedge clk or posedge rst) begin
		if (rst) begin
            valid2_r <= 0;
            for (n = 0; n < POINTS; n = n + 1) stage2_r[n] <= 0;
		end
		else if (valid1_r) begin
            valid2_r <= valid1_r;
            for (n = 0; n < POINTS/4; n = n + 1) begin
                a = n; b = n + POINTS/4;
                stage2_r[a] <= stage1_r[a] + stage1_r[b];
                stage2_r[b] <= stage1_r[a] - stage1_r[b];
            end

            for (n = 0; n < POINTS/4; n = n + 1) begin 
                c = POINTS/2 + n; d = POINTS/2 + n + POINTS/4;
                stage2_r[c] <= stage1_r[c] + stage1_i[d];
                diff_imag = stage1_r[c] + stage1_i[d];
                rot_r = diff_imag[7:0];
                stage2_r[d] <= rot_r;
            end
		end
	end
	
	 always @(posedge clk or posedge rst) begin
		if (rst) begin
            valid2_i <= 0;
            for (n = 0; n < POINTS; n = n + 1) stage2_i[n] <= 0;
		end
		else if (valid1_i) begin
            valid2_i <= valid1_i;
            for (n = 0; n < POINTS/4; n = n + 1) begin
                a = n; b = n + POINTS/4;
                stage2_i[a] <= stage1_i[a] + stage1_i[b];
                stage2_i[b] <= stage1_i[a] - stage1_i[b];
            end

            for (n = 0; n < POINTS/4; n = n + 1) begin 
                c = POINTS/2 + n; d = POINTS/2 + n + POINTS/4;
                stage2_i[c] <= stage1_i[c] - stage1_r[d];
                diff_real = stage1_i[c] - stage1_r[d];
                rot_i = -diff_real[7:0];
                stage2_i[d] <= rot_i;
            end
		end
	end

    // Stage 3 : 8-point butterflies
    always @(posedge clk or posedge rst) begin
        if(rst) begin
            valid3_r <= 0;
            for(n=0;n<POINTS;n=n+1) stage3_r[n]<=0;
        end else if(valid2_r) begin
            valid3_r <= valid2_r;
            for(n=0;n<POINTS/2;n=n+1) begin
                twiddle_idx = n;
                temp_r = stage2_r[2*n+1]*W_real[twiddle_idx] - stage2_i[2*n+1]*W_imag[twiddle_idx];
                mul_r = temp_r >>> SCALE;
                stage3_r[2*n]   <= stage2_r[2*n] + mul_r;
                stage3_r[2*n+1] <= stage2_r[2*n] - mul_r;
            end
        end
    end

    always @(posedge clk or posedge rst) begin
        if(rst) begin
            valid3_i <= 0;
            for(n=0;n<POINTS;n=n+1) stage3_i[n]<=0;
        end else if(valid2_i) begin
            valid3_i <= valid2_i;
            for(n=0;n<POINTS/2;n=n+1) begin
                twiddle_idx = n;
                temp_i = stage2_r[2*n+1]*W_imag[twiddle_idx] + stage2_i[2*n+1]*W_real[twiddle_idx];
                mul_i = temp_i >>> SCALE;
                stage3_i[2*n]   <= stage2_i[2*n] + mul_i;
                stage3_i[2*n+1] <= stage2_i[2*n] - mul_i;
            end
        end
    end

    // Output mapping
    always @(posedge clk) begin
        for(n=0;n<POINTS;n=n+1)
            Yr_flat[8*(n+1)-1 -: 8] <= stage3_r[n];
    end
	 
    always @(posedge clk) begin
        for(n=0;n<POINTS;n=n+1)
            Yi_flat[8*(n+1)-1 -: 8] <= stage3_i[n];
    end

    always @(posedge clk or posedge rst) begin
        if(rst) valid_out <= 0;
        else valid_out <= valid3_r;
    end
	 
endmodule
