class pf_description:
    def __init__( self,
                  path,
                  nonwater_tpr,
                  echo_trjconv,
                  initial_frame_trjconv,
                  final_frame_trjconv,
                  nonwater_xtc,
                  echo_rmsf,
                  initial_frame_rmsf,
                  final_frame_rmsf,
                  title_output):
        self.path = path
        self.nonwater_tpr = nonwater_tpr
        self.echo_trjconv = echo_trjconv
        self.initial_frame_trjconv = initial_frame_trjconv
        self.final_frame_trjconv = final_frame_trjconv
        self.nonwater_xtc = nonwater_xtc
        self.echo_rmsf = echo_rmsf
        self.initial_frame_rmsf = initial_frame_rmsf
        self.final_frame_rmsf = final_frame_rmsf
        self.title_output = title_output


    def get_path( self ):
        return self.path

    def get_nonwater_tpr( self ):
        return self.nonwater_tpr

    def get_echo_trjconv( self ):
        return self.echo_trjconv

    def get_initial_frame_trjconv( self ):
        return self.initial_frame_trjconv

    def get_final_frame_trjconv( self ):
        return self.final_frame_trjconv

    def get_nonwater_xtc(self):
        return self.nonwater_xtc

    def get_echo_rmsf(self):
        return self.echo_rmsf

    def get_initial_frame_rmsf(self):
        return self.initial_frame_rmsf

    def get_final_frame_rmsf(self):
        return self.final_frame_rmsf

    def get_title_output(self):
        return self.title_output

    def get_output_sufix_trjconv(self):

        output = str(self.get_initial_frame_trjconv()/1000) + "-" + str(self.get_final_frame_trjconv()/1000)

        if self.get_final_frame_trjconv() == self.get_final_frame_trjconv():
            output = str(self.get_initial_frame_trjconv()/1000)

        return output

    def get_output_sufix_rmsf( self ):

        output = str(self.get_initial_frame_rmsf()/1000) + "-" + str(self.get_final_frame_rmsf()/1000)

        if self.get_final_frame_rmsf() == self.get_final_frame_rmsf():
            output = str(self.get_initial_frame_rmsf()/1000)

        return output
