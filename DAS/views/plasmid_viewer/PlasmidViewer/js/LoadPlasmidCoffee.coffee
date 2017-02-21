#/**
# * Created by Justin on 2/13/17.
# */

'use strict'

#// Declare app level module which depends on views, and components

app = angular.module('myApp', ['angularplasmid'])


feature_loc = (feature, start, end, length) ->
  s = start(feature)
  e = end(feature)
  console.log(s, length, e)
  if s > length
    s = s - length
  if e > length
    e = e - length
  feature.plasmid_start = s
  feature.plasmid_end = e
  [s, e]

sanitize_features = (features, start, end, length) ->
  for f in features
    feature_loc(f, start, end, length)

get_feature_length = (a, length) ->
  x = a.plasmid_end - a.plasmid_start
  if x < 0
    x = x + length
  if x == NaN
    console.log('NNNN', a.plasmid_start, a.plasmid_end)
  x

sort_features = (features, length) ->
  # assign larget features first
  features.sort((a, b) ->
    get_feature_length(a, length) - get_feature_length(b, length)
  ).reverse()


assign_features = (features, length) ->

  feature_overlaps = (f1, f2, length) ->
    origin = false
    overlap = false
    s = f2.plasmid_start
    e = f2.plasmid_end
    if f1.shell == f2.shell
      if s > e
        origin = true

      if not origin
        if s <= f1.plasmid_start < e
          overlap = true
        if s < f1.plasmid_end <= e
          overlap = true
      else
        if s <= f1.plasmid_start <= length or 0 < f1.plasmid_start < e
          overlap = true
        if s < f1.plasmid_end <= length or 0 <= f1.plasmid_end <= e
          overlap = true
      if overlap
        console.log('\t' + f1.name + ' (' + f1.plasmid_start + '-' + f1.plasmid_end +  ') overlaps with ' + f2.name + ' (' + f2.plasmid_start + '-' + f2.plasmid_end + ') in shell ' + f1.shell)
    overlap

  max_feature_shell = 0
  assigned_features = []

  for feature in features
    if not feature.shell
      feature.shell = 1
      for af in assigned_features
        if feature_overlaps(feature, af) or feature_overlaps(af, feature)
          feature.shell = feature.shell + 1
      if feature.shell > max_feature_shell
        max_feature_shell = feature.shell
      console.log("Assigning " + feature.name + ' (' + feature.plasmid_start + '-' + feature.plasmid_end + ' to ' + feature.shell)
      assigned_features.push(feature)
  max_feature_shell



app.controller('LoadPlasmids', ($scope, $http) ->
#

  $http({method: 'POST', url: 'js/plasmiddata.json'}).success(
    (data) ->
      $scope.plasmid = data
      $scope.plasmid_size = 150
      $scope.feature_width = 10
      $scope.feature_padding = 3
      start = (x) -> x.start
      end = (x) -> x.stop

#
      sanitize_features($scope.plasmid.features, start, end, data.length)
      console.log("Sorting features")
      data.features = sort_features(data.features, data.length)
      for f in data.features
        console.log("L " + get_feature_length(f, data.length), f.plasmid_start, f.plasmid_end, data.length)
      $scope.max_feature_shell = assign_features($scope.plasmid.features, data.length)
#
      padding = 30
      $scope.contig_track_size = padding + $scope.plasmid_size + $scope.max_feature_shell*($scope.feature_width + $scope.feature_padding)
  )

  $http({method: 'POST', url: '../../alignment_viewer/data.json'}).success((data) ->
    $scope.contigs = data
    $scope.contig_padding = 10
    start = (x) -> x.query.start
    end = (x) -> x.query.end
    sanitize_features(data.contigs, start, end, $scope.plasmid.length)
    data.contigs = sort_features(data.features, data.length)
    $scope.max_contig_shell = assign_features(data.contigs, data.length)
  )

  $http({method: 'POST', url: '../../alignment_viewer/primer_data.json'}).success((data) -> $scope.primers = data)

  $http({method: 'POST', url: '../../alignment_viewer/bestassembly.json'}).success((data) -> $scope.assembly = data)
  $scope.color = () ->
    c = "rgba(" + (Math.ceil(Math.random() * 255)) + "," +
      (Math.ceil(Math.random() * 255)) + "," +
      (Math.ceil(Math.random() * 255)) + "," +
      0.9 + ")"
    c
)