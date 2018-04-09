import numpy as np

# load the pillar binary file, and wrap up into a class.
class PillarTrackerBinReader:
    fileversion1 = 589450  # file version control
    fileversion2 = 589451  # since v1.1.3 if grid restoration applied
    dt_int = np.dtype(np.int32).newbyteorder('>')
    dt_double = np.dtype(np.float64).newbyteorder('>')
    dt_int8 = np.dtype(np.int8).newbyteorder('>')

    def loadfile(self, fpath):
        fileversion1 = self.fileversion1
        fileversion2 = self.fileversion2
        dt_int = self.dt_int
        dt_double = self.dt_double
        dt_int8 = self.dt_int8

        class Meta:
            pass

        class ParaKernel:
            pass

        class ParaFFT:
            pass

        class Data:
            pass

        header = Meta()
        kernel = ParaKernel()
        fft = ParaFFT()
        data = Data()

        suc_load = False
        npillars = 0
        nframes = 0
        with open(fpath, "rb") as f:
            version = np.fromfile(f, dtype=dt_int, count=1)
            if version == fileversion1 or version == fileversion2:
                [npillars, nframes] = np.fromfile(f, dtype=dt_int, count=2)
                # nframes     = np.fromfile(f, dtype=dt_int, count=1)
                buffer = np.fromfile(f, dtype=dt_double, count=7)
                header.pixel_size = buffer[0]  # = np.fromfile(f, dtype=dt_double, count=1)
                header.diameter = buffer[1]  # = np.fromfile(f, dtype=dt_double, count=1)
                header.spacing = buffer[2]  # = np.fromfile(f, dtype=dt_double, count=1)
                header.oblique = buffer[3]  # = np.fromfile(f, dtype=dt_double, count=1)
                header.grid_angle = buffer[4]  # = np.fromfile(f, dtype=dt_double, count=1)
                header.sigma_PSF = buffer[5]  # = np.fromfile(f, dtype=dt_double, count=1)
                header.catch_radius = buffer[6]  # = np.fromfile(f, dtype=dt_double, count=1)
                [kernel_w, box_limit] = np.fromfile(f, dtype=dt_int, count=2)
                # kernel_w    = np.fromfile(f, dtype=dt_int, count=1)
                # box_limit   = np.fromfile(f, dtype=dt_int, count=1)
                header.kernel_w = kernel_w
                header.box_limit = box_limit

                buffer = np.fromfile(f, dtype=dt_int8, count=4)
                header.dark_pillar = buffer[0]  # np.fromfile(f, dtype=dt_int8, count=1)
                header.use_minimum_std = buffer[1]  # np.fromfile(f, dtype=dt_int8, count=1)
                header.use_metric_CG = buffer[2]  # np.fromfile(f, dtype=dt_int8, count=1)
                use_enhancer = buffer[3]  # np.fromfile(f, dtype=dt_int8, count=1)
                if version == fileversion2:
                    use_fft = np.fromfile(f, dtype=dt_int8, count=1)
                else:
                    use_fft = 0
                header.use_fft = use_fft
                header.use_enhancer = use_enhancer

                enhance_loaded_suc = True
                if use_enhancer > 0:
                    buffer = np.fromfile(f, dtype=dt_int8, count=3)
                    kernel.normalized_cc = buffer[0]  # np.fromfile(f, dtype=dt_int8, count=1)
                    kernel.use_gauss_psf = buffer[1]  # np.fromfile(f, dtype=dt_int8, count=1)
                    kernel.gauss_psf_dark = buffer[2]  # np.fromfile(f, dtype=dt_int8, count=1)
                    kernel.gauss_psf_radius = np.fromfile(f, dtype=dt_int, count=1)
                    kernel.gauss_psf_sigma = np.fromfile(f, dtype=dt_double, count=1)
                    [psf_w, psf_h] = np.fromfile(f, dtype=dt_int, count=1)
                    kernel.psf_w = psf_w
                    kernel.psf_h = psf_h
                    # psf_w = np.fromfile(f, dtype=dt_int, count=1)
                    # psf_h = np.fromfile(f, dtype=dt_int, count=1)
                    if psf_w > 0 and psf_h > 0:
                        PSF = np.fromfile(f, dtype=dt_double, count=psf_w * psf_h)
                        kernel.PSF = PSF.reshape(psf_h, psf_w)
                    else:
                        enhance_loaded_suc = False
                    header.kernel = kernel
                header.enhance_loaded_suc = enhance_loaded_suc

                if use_fft > 0:
                    fft_para = np.fromfile(f, dtype=dt_int, count=5)
                    fft.off_center_radius = fft_para[0]
                    fft.center_radius = fft_para[1]
                    fft.start_off_center_radius = fft_para[2]
                    fft.end_off_center_center_radius = fft_para[3]
                    num_points = fft_para[4]
                    if num_points > 0:
                        points = np.fromfile(f, dtype=dt_int, count=2 * num_points)
                        fft.points = points.reshape(2, num_points)
                    fft.num_points = num_points
                    header.fft = fft

                if enhance_loaded_suc:
                    startXY = np.fromfile(f, dtype=dt_int, count=2 * npillars)
                    header.startXY = startXY.reshape(2, npillars)
                    suc_load = True
            elif version % 2 == 0:
                npillars = version / 2
                nframes = np.fromfile(f, dtype=dt_int, count=1)
                suc_load = True
            # check whether the file header is loaded successfully.
            if suc_load and npillars > 0 and nframes > 0:
                ncols = npillars * 2
                nnsize = ncols * nframes
                A = np.fromfile(f, dtype=dt_double, count=nnsize)
                B = np.fromfile(f, dtype=dt_double, count=nnsize)
                A = A.reshape(ncols, nframes)
                B = B.reshape(ncols, nframes)
                end = ncols + 1;
                if version == fileversion2:
                    data.DX = B[0:end:2, :].transpose()
                    data.DY = B[1:end:2, :].transpose()
                    B = np.fromfile(f, dtype=dt_double, count=nnsize)
                    B = B.reshape(ncols, nframes)
                data.trackX = A[0:end:2, :].transpose()
                data.trackY = A[1:end:2, :].transpose()

                if B.size > 0:
                    data.drift_correctX = B[0:end:2, :].transpose()
                    data.drift_correctY = B[1:end:2, :].transpose()
                    driftXY = np.fromfile(f, dtype=dt_double, count=2 * nframes)
                    data.driftXY = driftXY.reshape(2, nframes)
                    data.reference_flags = np.fromfile(f, dtype=dt_int8, count=npillars)
            else:
                suc_load = False
                print('loading file failed.')

        header.suc_load = suc_load
        header.npillars = npillars
        header.nframes = nframes
        header.version = version
        return [header, data]
