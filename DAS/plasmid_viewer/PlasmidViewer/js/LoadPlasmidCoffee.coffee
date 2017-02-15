#/**
# * Created by Justin on 2/13/17.
# */

'use strict'

#// Declare app level module which depends on views, and components

app = angular.module('myApp', ['angularplasmid'])

app.controller('LoadPlasmids', ($scope, $http) ->

  $scope.size = 200

  $http({method: 'POST', url: 'js/plasmiddata.json'}).success(
    (data) ->

      feature_overlaps = (f1, f2) ->
        console.log("feature overlap")
        overlap = false
        if f2.shell == f1.shell
          s = f2.start
          e = f2.stop
          origin = false
          if s > e #then this feature goes over the origin
            origin = true
            s = f2.stop
            e = f2.start

          if s <= f1.start <= e
            if origin
              overlap = false
            else
              overlap = true

          if s <= f1.stop <= e
            if origin
              overlap = false
            else
              overlap = true
          if overlap
            console.log(s, e, f1.start, f1.stop)
        overlap

      alert data.length
      for plasmid in data
        assigned_features = []
        plasmid.features = plasmid.features.sort((a, b) ->
          x = Math.abs(a.start - a.stop)
          y = Math.abs(b.start - b.stop)
          x - y
        )
        plasmid.features.reverse()
        console.log([Math.abs(a.start - a.stop) for a in plasmid.features])

        for feature in plasmid.features
          if not feature.shell
            feature.shell = 1
            for af in assigned_features
                if feature_overlaps(feature, af)
                  feature.shell = feature.shell + 1
            assigned_features.push(feature)
#
      $scope.plasmids = data
  )

  $http({method: 'POST', url: '../../alignment_viewer/data.json'}).success((data) -> $scope.contigs = data)

  $http({method: 'POST', url: '../../alignment_viewer/primer_data.json'}).success((data) -> $scope.primers = data)

  $http({method: 'POST', url: '../../alignment_viewer/bestassembly.json'}).success((data) -> $scope.assembly = data)
  $scope.color = () ->
    c = "rgba(" + (Math.ceil(Math.random() * 255)) + "," +
      (Math.ceil(Math.random() * 255)) + "," +
      (Math.ceil(Math.random() * 255)) + "," +
      0.9 + ")"
    c
)