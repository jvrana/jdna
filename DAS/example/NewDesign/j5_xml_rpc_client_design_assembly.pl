#!/usr/bin/perl
# j5_xml_rpc_client_design_assembly.pl
# Nathan Hillson 01-SEP-15

# *** Copyright Notice ***
#
# Copyright (c) 2010-2015, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy). All rights
# reserved.
#
# NOTICE. This software was developed under funding from the U.S. Department of Energy.
# As such, the U.S. Government has been granted for itself and others acting on its
# behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, prepare derivative works, and perform publicly and display publicly.
# Beginning five (5) years after the date permission to assert copyright is obtained
# from the U.S. Department of Energy, and subject to any subsequent five (5) year
# renewals, the U.S. Government is granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce,
# prepare derivative works, distribute copies to the public, perform publicly and
# display publicly, and to permit others to do so. 
#
# This program is free software; you can redistribute it and/or modify
# it under the same terms as the Perl "Artistic License"
# (http://dev.perl.org/licenses/artistic.html).
#
# *** End of Copyright Notice ***

use XML::RPC;
use MIME::Base64;
require LWP::UserAgent;

my $assembly_method;
my $encoded_j5_parameters_file;
my $encoded_sequences_list_file;
my $encoded_zipped_sequences_file;
my $encoded_parts_list_file;
my $encoded_target_part_order_list_file;
my $encoded_eugene_rules_list_file;
my $encoded_master_plasmids_file;
my $encoded_master_oligos_file;
my $encoded_master_direct_syntheses_file;
my $input_file_handle;
my $j5_session_id;
my $master_direct_syntheses_list_filename;
my $master_oligos_list_filename;
my $master_plasmids_list_filename;
my $output_file_handle;
my $output_filename;
my $result;
my $reuse_j5_parameters_file;
my $reuse_sequences_list_file;
my $reuse_zipped_sequences_file;
my $reuse_parts_list_file;
my $reuse_target_part_order_list_file;
my $reuse_eugene_rules_list_file;
my $reuse_master_plasmids_file;
my $reuse_master_oligos_file;
my $reuse_master_direct_syntheses_file;
my $ua;
my $username;
my $xmlrpc;

if (!(@ARGV == 11 || (@ARGV == 9 && $ARGV[8]=~/Mock/))) {
    die "Usage:\n$0 j5_session_id TRUE_or_j5ParametersFileName TRUE_or_SequencesListFileName TRUE_or_ZippedSequencesFileName TRUE_or_PartsListFileName TRUE_or_TargetPartOrderListFileName TRUE_or_EugeneRulesListFileName TRUE_or_MasterPlasmidsListFileName TRUE_or_MasterOligosListFileName TRUE_or_MasterDirectSynthesesListFileName SLIC/Gibson/CPEC_or_GoldenGate_or_CombinatorialSLICGibsonCPEC_or_CombinatorialGoldenGate\nor\n$0 j5_session_id TRUE_or_j5ParametersFileName TRUE_or_SequencesListFileName TRUE_or_ZippedSequencesFileName TRUE_or_PartsListFileName TRUE_or_TargetPartOrderListFileName TRUE_or_EugeneRulesListFileName TRUE_or_MasterPlasmidsListFileName Mock_or_CombinatorialMock\n";
}

$j5_session_id = $ARGV[0];

# setting the user agent timeout to a longer period of time is required to make sure that j5 has sufficient time to complete an assembly
$ua = LWP::UserAgent->new();
$ua->timeout( 86400 ); # wait up until 24 hours for j5 to complete assembly, then bail
$xmlrpc = XML::RPC->new('https://j5.jbei.org/bin/j5_xml_rpc.pl',lwp_useragent => $ua);

# call the j5 xml-rpc service to see if last updated files are available
$result = $xmlrpc->call( 'GetLastUpdatedUserFiles', { j5_session_id => $j5_session_id } );

if (defined $result->{"faultString"}) {
  die "Error (GetLastUpdatedUserFiles): faultCode: ".$result->{"faultCode"}." faultString: ".$result->{"faultString"}."\n";
}

$username = $result->{"username"};

# initialize the parameters for calling the DesignAssembly method
# j5 parameters file
if ($ARGV[1]=~/^TRUE$/i) {
  if ($result->{"j5_parameters_file_exists"}=~/^TRUE$/i) {
    $reuse_j5_parameters_file="TRUE";
  }
  else {
    die "Error (DesignAssembly): j5 parameters file not found on server for user $username";
  }
}
else {
  $reuse_j5_parameters_file="FALSE";

  # encode the j5_parameters_file
  open($input_file_handle, "$ARGV[1]") or die "$!";
  while (read($input_file_handle, $input_buffer, 60*57)) {
    $encoded_j5_parameters_file.=encode_base64($input_buffer);
  }
  chomp $encoded_j5_parameters_file;
  close($input_file_handle);
}

# sequences list file
if ($ARGV[2]=~/^TRUE$/i) {
  if ($result->{"sequences_list_file_exists"}=~/^TRUE$/i) {
    $reuse_sequences_list_file="TRUE";
  }
  else {
    die "Error (DesignAssembly): sequences list file not found on server for user $username";
  }
}
else {
  $reuse_sequences_list_file="FALSE";

  # encode the sequences_list_file
  open($input_file_handle, "$ARGV[2]") or die "$!";
  while (read($input_file_handle, $input_buffer, 60*57)) {
    $encoded_sequences_list_file.=encode_base64($input_buffer);
  }
  chomp $encoded_sequences_list_file;
  close($input_file_handle);
}

# zipped sequences list file
if ($ARGV[3]=~/^TRUE$/i) {
  if ($result->{"zipped_sequences_file_exists"}=~/^TRUE$/i) {
    $reuse_zipped_sequences_file="TRUE";
  }
  else {
    die "Error (DesignAssembly): zipped sequences file not found on server for user $username";
  }
}
else {
  $reuse_zipped_sequences_file="FALSE";

  # encode the zipped_sequences_file
  open($input_file_handle, "$ARGV[3]") or die "$!";
  while (read($input_file_handle, $input_buffer, 60*57)) {
    $encoded_zipped_sequences_file.=encode_base64($input_buffer);
  }
  chomp $encoded_zipped_sequences_file;
  close($input_file_handle);
}

# parts list file
if ($ARGV[4]=~/^TRUE$/i) {
  if ($result->{"parts_list_file_exists"}=~/^TRUE$/i) {
    $reuse_parts_list_file="TRUE";
  }
  else {
    die "Error (DesignAssembly): parts list file not found on server for user $username";
  }
}
else {
  $reuse_parts_list_file="FALSE";

  # encode the parts_list_file
  open($input_file_handle, "$ARGV[4]") or die "$!";
  while (read($input_file_handle, $input_buffer, 60*57)) {
    $encoded_parts_list_file.=encode_base64($input_buffer);
  }
  chomp $encoded_parts_list_file;
  close($input_file_handle);
}

# target part order list file
if ($ARGV[5]=~/^TRUE$/i) {
  if ($result->{"target_part_order_list_file_exists"}=~/^TRUE$/i) {
    $reuse_target_part_order_list_file="TRUE";
  }
  else {
    die "Error (DesignAssembly): target part order list file not found on server for user $username";
  }
}
else {
  $reuse_target_part_order_list_file="FALSE";

  # encode the target_part_order_list_file
  open($input_file_handle, "$ARGV[5]") or die "$!";
  while (read($input_file_handle, $input_buffer, 60*57)) {
    $encoded_target_part_order_list_file.=encode_base64($input_buffer);
  }
  chomp $encoded_target_part_order_list_file;
  close($input_file_handle);
}

# eugene rules list file
if ($ARGV[6]=~/^TRUE$/i) {
  if ($result->{"eugene_rules_list_file_exists"}=~/^TRUE$/i) {
    $reuse_eugene_rules_list_file="TRUE";
  }
  else {
    die "Error (DesignAssembly): eugene rules list file not found on server for user $username";
  }
}
else {
  $reuse_eugene_rules_list_file="FALSE";

  # encode the eugene_rules_list_file
  open($input_file_handle, "$ARGV[6]") or die "$!";
  while (read($input_file_handle, $input_buffer, 60*57)) {
    $encoded_eugene_rules_list_file.=encode_base64($input_buffer);
  }
  chomp $encoded_eugene_rules_list_file;
  close($input_file_handle);
}

# master plasmids file
if ($ARGV[7]=~/^TRUE$/i) {
  if ($result->{"master_plasmids_file_exists"}=~/^TRUE$/i) {
    $reuse_master_plasmids_file="TRUE";
  }
  else {
    die "Error (DesignAssembly): master plasmids file not found on server for user $username";
  }
}
else {
  $reuse_master_plasmids_file="FALSE";
  $master_plasmids_list_filename = $ARGV[7];

  # encode the master_plasmids_file
  open($input_file_handle, "$ARGV[7]") or die "$!";
  while (read($input_file_handle, $input_buffer, 60*57)) {
    $encoded_master_plasmids_file.=encode_base64($input_buffer);
  }
  chomp $encoded_master_plasmids_file;
  close($input_file_handle);
}

if (@ARGV==11) {
  # master oligos file
  if ($ARGV[8]=~/^TRUE$/i) {
    if ($result->{"master_oligos_file_exists"}=~/^TRUE$/i) {
      $reuse_master_oligos_file="TRUE";
    }
    else {
      die "Error (DesignAssembly): master oligos file not found on server for user $username";
    }
  }
  else {
    $reuse_master_oligos_file="FALSE";
    $master_oligos_list_filename = $ARGV[8];
    
    # encode the master_oligos_list_file
    open($input_file_handle, "$ARGV[8]") or die "$!";
    while (read($input_file_handle, $input_buffer, 60*57)) {
      $encoded_master_oligos_file.=encode_base64($input_buffer);
    }
    chomp $encoded_master_oligos_file;
    close($input_file_handle);
  }
  
  # master direct syntheses file
  if ($ARGV[9]=~/^TRUE$/i) {
    if ($result->{"master_direct_syntheses_file_exists"}=~/^TRUE$/i) {
      $reuse_master_direct_syntheses_file="TRUE";
    }
    else {
      die "Error (DesignAssembly): master direct syntheses file not found on server for user $username";
    }
  }
  else {
    $reuse_master_direct_syntheses_file="FALSE";
    $master_direct_syntheses_list_filename = $ARGV[9];
    
    # encode the master_direct_syntheses_list_file
    open($input_file_handle, "$ARGV[9]") or die "$!";
    while (read($input_file_handle, $input_buffer, 60*57)) {
      $encoded_master_direct_syntheses_file.=encode_base64($input_buffer);
    }
    chomp $encoded_direct_syntheses_file;
    close($input_file_handle);
  }

  $assembly_method = $ARGV[10];

  # call the j5 xml-rpc service to design the assembly
  $result = $xmlrpc->call( 'DesignAssembly', { 
					      j5_session_id => $j5_session_id,
					      reuse_j5_parameters_file => $reuse_j5_parameters_file,
					      encoded_j5_parameters_file => $encoded_j5_parameters_file,
					      reuse_sequences_list_file => $reuse_sequences_list_file,
					      encoded_sequences_list_file => $encoded_sequences_list_file,
					      reuse_zipped_sequences_file => $reuse_zipped_sequences_file,
					      encoded_zipped_sequences_file => $encoded_zipped_sequences_file,
					      reuse_parts_list_file => $reuse_parts_list_file,
					      encoded_parts_list_file => $encoded_parts_list_file,
					      reuse_target_part_order_list_file => $reuse_target_part_order_list_file,
					      encoded_target_part_order_list_file => $encoded_target_part_order_list_file,
					      reuse_eugene_rules_list_file => $reuse_eugene_rules_list_file,
					      encoded_eugene_rules_list_file => $encoded_eugene_rules_list_file,
					      reuse_master_plasmids_file => $reuse_master_plasmids_file,
					      encoded_master_plasmids_file => $encoded_master_plasmids_file,
					      master_plasmids_list_filename => $master_plasmids_list_filename,
					      reuse_master_oligos_file => $reuse_master_oligos_file,
					      encoded_master_oligos_file => $encoded_master_oligos_file,
					      master_oligos_list_filename => $master_oligos_list_filename,
					      reuse_master_direct_syntheses_file => $reuse_master_direct_syntheses_file,
					      encoded_master_direct_syntheses_file => $encoded_master_direct_syntheses_file,
					      master_direct_syntheses_list_filename => $master_direct_syntheses_list_filename,
					      assembly_method => $assembly_method,
					     }
			 );
}
else {
  $assembly_method = $ARGV[8];

  # call the j5 xml-rpc service to design the assembly
  $result = $xmlrpc->call( 'DesignAssembly', { 
					      j5_session_id => $j5_session_id,
					      reuse_j5_parameters_file => $reuse_j5_parameters_file,
					      encoded_j5_parameters_file => $encoded_j5_parameters_file,
					      reuse_sequences_list_file => $reuse_sequences_list_file,
					      encoded_sequences_list_file => $encoded_sequences_list_file,
					      reuse_zipped_sequences_file => $reuse_zipped_sequences_file,
					      encoded_zipped_sequences_file => $encoded_zipped_sequences_file,
					      reuse_parts_list_file => $reuse_parts_list_file,
					      encoded_parts_list_file => $encoded_parts_list_file,
					      reuse_target_part_order_list_file => $reuse_target_part_order_list_file,
					      encoded_target_part_order_list_file => $encoded_target_part_order_list_file,
					      reuse_eugene_rules_list_file => $reuse_eugene_rules_list_file,
					      encoded_eugene_rules_list_file => $encoded_eugene_rules_list_file,
					      reuse_master_plasmids_file => $reuse_master_plasmids_file,
					      encoded_master_plasmids_file => $encoded_master_plasmids_file,
					      master_plasmids_list_filename => $master_plasmids_list_filename,
					      assembly_method => $assembly_method,
					     }
		       );
}

if (defined $result->{"faultString"}) {
  print "Error: faultCode: ",$result->{"faultCode"}," faultString: ",$result->{"faultString"},"\n";
}

else {

  if (defined $result->{"error_message"}) {
    print "Error: ",$result->{"error_message"},"\n";
  }

  $output_filename = $result->{"output_filename"};
  open($output_file_handle, ">$output_filename") or die "$!";
  print $output_file_handle decode_base64($result->{"encoded_output_file"});
  close($output_file_handle);
}

#print $xmlrpc->xml_in();
